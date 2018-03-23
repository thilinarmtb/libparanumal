#include "multiFluidIns2D.h"

multiFluidIns_t *multiFluidInsLevelSetSetup2D(mesh2D *mesh, char * options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", (rank)%3);
  //  sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
  //printf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 0");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
  //sprintf(deviceConfig, "mode = Serial");  

  multiFluidIns_t *levelSet = (multiFluidIns_t*) calloc(1, sizeof(multiFluidIns_t));


  levelSet->Nfields  = 1; // Each Velocity Field
  levelSet->index    = 0; // no history for velocity
  
  mesh->Nfields = levelSet->Nfields;

  levelSet->mesh = mesh;

  int Ntotal = mesh->Np*mesh->Nelements;
  levelSet->Nblock = (Ntotal+blockSize-1)/blockSize;

  int Nstages = 3;
  // compute samples of q at interpolation nodes
  levelSet->U      = (dfloat*) calloc((mesh->totalHaloPairs + mesh->Nelements)*mesh->Np,sizeof(dfloat));
  levelSet->V      = (dfloat*) calloc((mesh->totalHaloPairs + mesh->Nelements)*mesh->Np,sizeof(dfloat));
  levelSet->Phi    = (dfloat*) calloc((mesh->totalHaloPairs + mesh->Nelements)*mesh->Np,sizeof(dfloat));
  levelSet->rhsPhi = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  // Fill up required fileds
  levelSet->finalTime = 1.0;

  // Initial Conditions, Flow Properties
  if (rank==0) 
    printf("Starting initial conditions for LevelSet2D\n");

  dfloat ux   = 0.0  ;
  dfloat uy   = 0.0  ;
  dfloat phi  = 0.0  ;
  
  // Define total DOF per field for INS i.e. (Nelelemts + Nelements_halo)*Np
  levelSet->NtotalDofs = (mesh->totalHaloPairs+mesh->Nelements)*mesh->Np ;
  levelSet->NDofs      = mesh->Nelements*mesh->Np;
  // Initialize
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      const int id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      #if 1  
        dfloat xc = 0.50; 
        dfloat yc = 0.75;
        dfloat r  = 0.15;  
        //   
        levelSet->U[id]   = -2.*M_PI*(y-0.5);
        levelSet->V[id]   =  2.*M_PI*(x-0.5);
        levelSet->Phi[id] = sqrt( (x-xc)*(x-xc)+(y-yc)*(y-yc)) - r;
      #endif
    }
  }

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }

  // Find Maximum Velocity
  dfloat umax = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      const int id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = levelSet->U[id];
      dfloat uyn = levelSet->V[id];

      //Squared maximum velocity
      dfloat numax = uxn*uxn + uyn*uyn;
      umax = mymax(umax, numax);
    }
  }
  // Maximum Velocity
  umax = sqrt(umax);

 
  dfloat cfl = 0.5; 
 
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  dfloat dt     = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  if (rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  cfl);
    printf("dt = %g\n",   dt);
  }

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(levelSet->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  levelSet->NtimeSteps = levelSet->finalTime/levelSet->dt;
  levelSet->dt         = levelSet->finalTime/levelSet->NtimeSteps;

  levelSet->errorStep = 100;

  if (rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", levelSet->NtimeSteps,levelSet->errorStep, levelSet->dt);

  occa::kernelInfo kernelInfo;
  meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo.addDefine("p_maxNodesVolumeCub", maxNodesVolumeCub);

  int cubNblockV = mymax(1,256/maxNodesVolumeCub); 
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxNodesSurfaceCub",maxNodesSurfaceCub);

  int cubNblockS = mymax(1,256/maxNodesSurfaceCub); // works for CUDA
  //

  if (rank==0) printf("p_cubNblockV=%d\n", cubNblockV);
  if (rank==0) printf("p_cubNblockS=%d\n", cubNblockS);
  kernelInfo.addDefine("p_cubNblockV",cubNblockV);
  kernelInfo.addDefine("p_cubNblockS",cubNblockS);



  if (rank==0) printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  if (rank==0) printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

  if (rank==0) printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);

  // ADD-DEFINES
  // kernelInfo.addDefine("p_Lambda2", 0.5f);
  // kernelInfo.addDefine("p_NTfields", levelSet->NTfields);
  // kernelInfo.addDefine("p_NVfields", levelSet->NVfields);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfaces*mesh->Nfp);
  
 
  levelSet->o_U       = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), levelSet->U);
  levelSet->o_V       = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), levelSet->V);
  levelSet->o_Phi     = mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), levelSet->Phi);
  levelSet->o_rhsPhi  = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), levelSet->rhsPhi);
  levelSet->o_resPhi  = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), levelSet->rhsPhi);

  if(mesh->totalHaloPairs){
    levelSet->o_phiHaloBuffer = mesh->device.malloc(mesh->totalHaloPairs*mesh->Np *sizeof(dfloat));
  }

  levelSet->mesh = mesh;




  // ===========================================================================
  levelSet->levelSetHaloExtractKernel=
      mesh->device.buildKernelFromSource(DHOLMES "/okl/multiFluidInsHaloExchange.okl",
        "multiFluidInsLevelSetHaloExtract",
          kernelInfo);

  levelSet->levelSetHaloScatterKernel=
      mesh->device.buildKernelFromSource(DHOLMES "/okl/multiFluidInsHaloExchange.okl",
        "multiFluidInsLevelSetHaloScatter",
          kernelInfo); 

   levelSet->levelSetTransportVolumeKernel = 
      mesh->device.buildKernelFromSource(DHOLMES "/okl/multiFluidInsLevelSetTransport.okl",
        "multiFluidInsLevelSetTransportCubatureVolume2D",
          kernelInfo);  

    levelSet->levelSetTransportSurfaceKernel = 
          mesh->device.buildKernelFromSource(DHOLMES "/okl/multiFluidInsLevelSetTransport.okl",
            "multiFluidInsLevelSetTransportCubatureSurface2D",
              kernelInfo); 

    levelSet->levelSetTransportUpdateKernel = 
          mesh->device.buildKernelFromSource(DHOLMES "/okl/multiFluidInsLevelSetTransport.okl",
            "multiFluidInsLevelSetTransportLSERKUpdate2D",
              kernelInfo);  

  // ins->scaledAddKernel =
  //     mesh->device.buildKernelFromSource(DHOLMES "/okl/scaledAdd.okl",
  //          "scaledAddwOffset",
  //          kernelInfo);
 
  // ins->totalHaloExtractKernel=
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //            "insTotalHaloExtract2D",
  //            kernelInfo);

  // ins->totalHaloScatterKernel=
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //            "insTotalHaloScatter2D",
  //            kernelInfo);
  // ins->velocityHaloExtractKernel=
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //           "insVelocityHaloExtract2D",
  //            kernelInfo);

  // ins->velocityHaloScatterKernel=
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //            "insVelocityHaloScatter2D",
  //            kernelInfo);

  // ins->pressureHaloExtractKernel=
  //   mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //            "insPressureHaloExtract",
  //            kernelInfo);

  // ins->pressureHaloScatterKernel=
  //  mesh->device.buildKernelFromSource(DHOLMES "/okl/insHaloExchange.okl",
  //            "insPressureHaloScatter",
  //            kernelInfo);  

  return levelSet;
}







