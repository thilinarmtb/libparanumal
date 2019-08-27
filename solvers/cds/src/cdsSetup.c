/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "cds.h"
#include "omp.h"
#include <unistd.h>

cds_t *cdsSetup(mesh_t *mesh, setupAide options){

  cds_t *cds = new cds_t();
  cds->mesh = mesh;
  cds->options = options;
  cds->kernelInfo = new occa::properties(); 

  options.getArgs("MESH DIMENSION", cds->dim);
  options.getArgs("ELEMENT TYPE", cds->elementType);

  cds->NVfields = (cds->dim ==3) ? 3:2; // Velocity Fields  
  cds->NSfields = 1 ;  //  Total Number of Scalar  Fields
  mesh->Nfields = 1; 
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    cds->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    cds->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    cds->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));
    cds->extC    = (dfloat*) calloc(3, sizeof(dfloat));
  }else if(!options.compareArgs("TIME INTEGRATOR", "EXTBDF")){
    printf("Currently only BDF time stepping is implemented\n");
    exit(-1);
  }
      
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    cds->Nstages = 1;
    cds->temporalOrder = 1;
    cds->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    cds->Nstages = 2;
    cds->temporalOrder = 2;
    cds->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    cds->Nstages = 3;
    cds->temporalOrder = 3;
    cds->g0 = 11.f/6.f;
  }


  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  cds->Ntotal      = Ntotal;
  cds->vOffset     = Ntotal; // keep it different for now.....
  cds->sOffset     = Ntotal;
  cds->Nblock      = (Nlocal+blockSize-1)/blockSize;

  // Solution storage at interpolation nodes
  cds->U     = (dfloat*) calloc(cds->NVfields*(cds->Nstages+0)*Ntotal,sizeof(dfloat));
  cds->S     = (dfloat*) calloc(cds->NSfields*(cds->Nstages+0)*Ntotal,sizeof(dfloat));
  cds->NS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*Ntotal,sizeof(dfloat));
  cds->rkS   = (dfloat*) calloc(cds->NSfields                 *Ntotal,sizeof(dfloat));
  cds->rhsS  = (dfloat*) calloc(cds->NSfields                 *Ntotal,sizeof(dfloat));
  
  cds->Nsubsteps = 0;
   
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",cds->Nsubsteps);


  if(cds->Nsubsteps){
    cds->Sd      = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
    cds->resS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
    cds->rhsSd   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));        
    cds->Ue      = (dfloat*) calloc(cds->NVfields*Ntotal,sizeof(dfloat));

 
    // Prepare RK stages for Subcycling Part
    
    int Sorder = 4; // Defaulting to LSERK 4(5) 
    
    options.getArgs("SUBCYCLING TIME ORDER", Sorder);

    if(Sorder==2){
      cds->SNrk     = 2; 
      dfloat rka[2] = {0.0,     1.0 };
      dfloat rkb[2] = {0.5,     0.5 };
      dfloat rkc[2] = {0.0,     1.0 };
      //
      cds->Srka = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkb = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkc = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      //
      memcpy(cds->Srka, rka, cds->SNrk*sizeof(dfloat));
      memcpy(cds->Srkb, rkb, cds->SNrk*sizeof(dfloat));
      memcpy(cds->Srkc, rkc, cds->SNrk*sizeof(dfloat));
    }else if(Sorder ==3){
      // Using Williamson 3rd order scheme converted to low storage since the better truncation 
      cds->SNrk     = 3; 
      dfloat rka[3] = {0.0,     -5.0/9.0,  -153.0/128.0};
      dfloat rkb[3] = {1.0/3.0, 15.0/16.0,    8.0/15.0 };
      dfloat rkc[3] = {0.0,      1.0/3.0,     3.0/4.0  };
      //
      cds->Srka = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkb = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkc = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      //
      memcpy(cds->Srka, rka, cds->SNrk*sizeof(dfloat));
      memcpy(cds->Srkb, rkb, cds->SNrk*sizeof(dfloat));
      memcpy(cds->Srkc, rkc, cds->SNrk*sizeof(dfloat));
    }else{
      cds->SNrk     = 5; 
      cds->Srka = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkb = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      cds->Srkc = (dfloat*) calloc(cds->SNrk, sizeof(dfloat));
      // Assumes initialized in mesh, can be moved here
      for(int rk=0; rk<cds->SNrk; rk++){
      	cds->Srka[rk] = mesh->rka[rk]; 
      	cds->Srkb[rk] = mesh->rkb[rk]; 
      	cds->Srkc[rk] = mesh->rkc[rk]; 
      }
    }
  }

  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", cds->ubar);
  options.getArgs("VBAR", cds->vbar);
  if (cds->dim==3)
    options.getArgs("WBAR", cds->wbar);
  options.getArgs("SBAR", cds->sbar);

  options.getArgs("CONDUCTIVITY", cds->k);
  options.getArgs("SPECIFIC HEAT", cds->cp);
  options.getArgs("DENSITY", cds->rho);
  options.getArgs("VISCOSITY", cds->nu);

  cds->alf    = cds->k/(cds->rho*cds->cp); // Thermal diffusivity
  cds->ialf   = 1.0/cds->alf;              // Inverse diff. 


  //Reynolds number
  cds->Re     = cds->ubar/cds->nu;         // Reynolds number
  cds->Pr     = cds->nu/cds->alf;          // Prandtl number
    
  occa::properties& kernelInfo = *cds->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(cds->dim==3){
    if(cds->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo); 
  } 
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  // ADD-DEFINES
  kernelInfo["defines/" "p_sbar"]= cds->sbar;
  kernelInfo["defines/" "p_ubar"]= cds->ubar;
  kernelInfo["defines/" "p_vbar"]= cds->vbar;
  kernelInfo["defines/" "p_wbar"]= cds->wbar;

  kernelInfo["defines/" "p_NSfields"]= cds->NSfields;
  kernelInfo["defines/" "p_NVfields"]= cds->NVfields;
  kernelInfo["defines/" "p_NTfields"]= (cds->NVfields+cds->NSfields);
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  cds->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  cds->Nsubsteps;


  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();
 
  cds->o_U = mesh->device.malloc(cds->NVfields*(cds->Nstages+0)*Ntotal*sizeof(dfloat), cds->U);
  cds->o_S = mesh->device.malloc(cds->NSfields*(cds->Nstages+0)*Ntotal*sizeof(dfloat), cds->S);

#if 0
  if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);
#endif
 
  for (int r=0;r<2;r++){
      if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
      if (cds->dim==2){ 
        cds->setFlowFieldKernel   =  mesh->device.buildKernel(DCDS "/okl/cdsSetFlowField2D.okl", "cdsSetFlowField2D", kernelInfo);
        cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField2D.okl", "cdsSetScalarField2D", kernelInfo);
      }else{
        cds->setFlowFieldKernel   =  mesh->device.buildKernel(DCDS "/okl/cdsSetFlowField3D.okl", "cdsSetFlowField3D", kernelInfo);  
        cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField3D.okl", "cdsSetScalarField3D", kernelInfo);  
      }
    }
    MPI_Barrier(mesh->comm);
  }

  cds->startTime =0.0;
  options.getArgs("START TIME", cds->startTime);

  cds->setFlowFieldKernel(mesh->Nelements,
                          cds->startTime,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          cds->vOffset,
                          cds->o_U);

  cds->setScalarFieldKernel(mesh->Nelements,
			    cds->startTime,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    cds->sOffset,
			    cds->o_S);

  cds->o_U.copyTo(cds->U);
  cds->o_S.copyTo(cds->S);


#if 0
  printf("AK: Writing initial condition file and exitiong ....\n")
  char fname[BUFSIZ];
    string outName;
    cds->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, cds->frame++);
    cdsPlotVTU(cds, fname);
    return cds; 
#endif




  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){

    if(cds->elementType==TRIANGLES || cds->elementType == TETRAHEDRA){
      for(int f=0;f<mesh->Nfaces;++f){
        dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
        dfloat sJ   = mesh->sgeo[sid + SJID];
        dfloat invJ = mesh->sgeo[sid + IJID];

        dfloat hest = 2./(sJ*invJ);

        hmin = mymin(hmin, hest);
        hmax = mymax(hmax, hest);
      }
    }else{
      for(int f=0;f<mesh->Nfaces;++f){
        for(int n=0; n<mesh->Nfp; n++){
	  dlong sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n);
	  dfloat sJ   = mesh->sgeo[sid + SJID];
	  dfloat invJ = mesh->sgeo[sid + IJID];

	  dfloat hest = 2./(sJ*invJ);

	  hmin = mymin(hmin, hest);
	  hmax = mymax(hmax, hest);
	}
      }
    }

    // dfloat maxMagVecLoc = 0;

    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = cds->U[id+0*cds->vOffset];
      dfloat uyn = cds->U[id+1*cds->vOffset];
      dfloat uzn = 0.0;
      if (cds->dim==3) uzn = cds->U[id+2*cds->vOffset];


      //Squared maximum velocity
      dfloat numax;
      if (cds->dim==2)
        numax = uxn*uxn + uyn*uyn;
      else 
        numax = uxn*uxn + uyn*uyn + uzn*uzn;

      umax = mymax(umax, numax);
    }
  }

  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity

  options.getArgs("CFL", cds->cfl);
  dfloat dt     = cds->cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
  
  cds->dtAdaptStep = 0; 
  options.getArgs("TSTEPS FOR TIME STEP ADAPT", cds->dtAdaptStep);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(cds->dti), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  // save initial time-step estimate 
  cds->dt = cds->dti; 

  options.getArgs("FINAL TIME", cds->finalTime);
  options.getArgs("START TIME", cds->startTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    cds->NtimeSteps = (cds->finalTime-cds->startTime)/cds->dt;

    if(cds->Nsubsteps){
      cds->dt         = cds->Nsubsteps*cds->dt;
      cds->NtimeSteps = (cds->finalTime-cds->startTime)/cds->dt;
      cds->dt         = (cds->finalTime-cds->startTime)/cds->NtimeSteps;
      cds->sdt        = cds->dt/cds->Nsubsteps;
    } else{
      cds->NtimeSteps = (cds->finalTime-cds->startTime)/cds->dt;
      cds->dt         = (cds->finalTime-cds->startTime)/cds->NtimeSteps;
    }
  }

  cds->dtMIN = 1E-2*cds->dt; //minumum allowed timestep

  if (mesh->rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  cds->cfl);
    printf("dt = %g\n",   cds->dt);
  }
 
  if (cds->Nsubsteps && mesh->rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", cds->dt, cds->sdt, cds->dt/cds->sdt);

  
  // Hold some inverses for kernels
  //cds->ialf  = 1.0/(cds->k/cds->rho*cds->cp); 
  cds->idt     = 1.0/cds->dt;
  cds->lambda  = cds->g0 / (cds->dt * cds->alf);

  kernelInfo["defines/" "p_alf"]     = cds->alf;
  // kernelInfo["defines/" "p_irhocp"]  = cds->ialf;
  // kernelInfo["defines/" "p_nu"]      = cds->nu;

  
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", cds->outputStep);
  if (mesh->rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", cds->NtimeSteps,cds->outputStep, cds->dt);

  cds->outputForceStep = 0;
  options.getArgs("TSTEPS FOR FORCE OUTPUT", cds->outputForceStep);

  //make option objects for elliptc solvers
  cds->options = options;
  cds->options.setArgs("KRYLOV SOLVER",        options.getArgs("KRYLOV SOLVER"));
  cds->options.setArgs("DISCRETIZATION",       options.getArgs("DISCRETIZATION"));
  cds->options.setArgs("BASIS",                options.getArgs("BASIS"));
  cds->options.setArgs("PRECONDITIONER",       options.getArgs("PRECONDITIONER"));
  cds->options.setArgs("MULTIGRID COARSENING", options.getArgs("MULTIGRID COARSENING"));
  cds->options.setArgs("MULTIGRID SMOOTHER",   options.getArgs("MULTIGRID SMOOTHER"));
  cds->options.setArgs("PARALMOND CYCLE",      options.getArgs("PARALMOND CYCLE"));
  cds->options.setArgs("PARALMOND SMOOTHER",   options.getArgs("PARALMOND SMOOTHER"));
  cds->options.setArgs("PARALMOND PARTITION",  options.getArgs("PARALMOND PARTITION"));

  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall bc = 2 -> inflow bc = 3 -> outflow
  // bc = 4 -> x-aligned slip, bc = 5 -> y-aligned slip, bc = 6 -> z-aligned slip

  int sBCType[7] = {0,1,1,2,1,1,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
 
  //Solver tolerances 
  cds->TOL = 1E-10;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  cds->solver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  cds->solver->mesh = mesh;
  cds->solver->options = cds->options;
  cds->solver->dim = cds->dim;
  cds->solver->elementType = cds->elementType;
  cds->solver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(cds->solver->BCType,sBCType,7*sizeof(int));
  ellipticSolveSetup(cds->solver, cds->lambda, kernelInfo); 

  //make node-wise boundary flags
  cds->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) cds->mapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          cds->mapB[fid+e*mesh->Np] = mymin(bc,cds->mapB[fid+e*mesh->Np]);
        }
      }
    }
  }

  ogsGatherScatter(cds->mapB, ogsInt, ogsMin, mesh->ogs);
  
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (cds->mapB[n] == 1E9) {
      cds->mapB[n] = 0.;
    }
  }
  cds->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), cds->mapB);

  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  if(options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    kernelInfo["defines/" "p_EXTBDF"]= 1;
  else
    kernelInfo["defines/" "p_EXTBDF"]= 0;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  


#if 1
  if (mesh->rank==0) {
    printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
    printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

    printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  }
#endif
  

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    cds->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    cds->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    cds->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    cds->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    cds->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    cds->o_prkA = cds->o_extbdfC;
    cds->o_prkB = cds->o_extbdfC;
  }else{
    printf("Only BDF is implemented\n");
    exit(EXIT_FAILURE);
  }

  // MEMORY ALLOCATION
  cds->o_rhsS  = mesh->device.malloc(cds->NSfields*                 Ntotal*sizeof(dfloat), cds->rhsS);
  cds->o_NS    = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*Ntotal*sizeof(dfloat), cds->NS);
  cds->o_rkS   = mesh->device.malloc(                               Ntotal*sizeof(dfloat), cds->rkS);  

  if(mesh->totalHaloPairs){//halo setup
    // Define new variable nodes per element (npe) to test thin halo; 
//    int npe = mesh->Np; // will be changed 
    int npe = mesh->Nfp; 
    dlong haloBytes   = mesh->totalHaloPairs*npe*(cds->NSfields + cds->NVfields)*sizeof(dfloat);
    dlong gatherBytes = (cds->NSfields+cds->NVfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
    cds->o_haloBuffer = mesh->device.malloc(haloBytes);
#if 0
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(gatherBytes, NULL);
    
    cds->sendBuffer    = (dfloat*) o_sendBuffer.getMappedPointer();
    cds->recvBuffer    = (dfloat*) o_recvBuffer.getMappedPointer();
    cds->haloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
#endif
    //occa::memory o_sendBuffer, o_recvBuffer,o_gatherTmpPinned;

    cds->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_sendBuffer, cds->h_sendBuffer);
    cds->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_recvBuffer, cds->h_recvBuffer);
    cds->haloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, gatherBytes, NULL, cds->o_gatherTmpPinned, cds->h_gatherTmpPinned); 
    cds->o_haloGatherTmp = mesh->device.malloc(gatherBytes,  cds->haloGatherTmp);

    // Halo exchange for more efficient subcycling 
    if(cds->Nsubsteps){
      dlong shaloBytes   = mesh->totalHaloPairs*npe*(cds->NSfields)*sizeof(dfloat);
      dlong sgatherBytes = (cds->NSfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
      cds->o_shaloBuffer = mesh->device.malloc(shaloBytes);

      //occa::memory o_ssendBuffer, o_srecvBuffer,o_sgatherTmpPinned;

      cds->ssendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL, cds->o_ssendBuffer, cds->h_ssendBuffer);
      cds->srecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL, cds->o_srecvBuffer, cds->h_srecvBuffer);
      cds->shaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, sgatherBytes, NULL, cds->o_sgatherTmpPinned, cds->h_sgatherTmpPinned);
      cds->o_shaloGatherTmp = mesh->device.malloc(sgatherBytes,  cds->shaloGatherTmp);
   }
  }

  // set kernel name suffix
  char *suffix;
  
  if(cds->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(cds->elementType==QUADRILATERALS){
    if(cds->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D"); 
  }
  if(cds->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(cds->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

   for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
      
      sprintf(fileName, DCDS "/okl/cdsHaloExchange.okl");

      
      // sprintf(kernelName, "cdsHaloExtract");
      // cds->haloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
      
      // sprintf(kernelName, "cdsHaloScatter");
      // cds->haloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // if(cds->Nsubsteps){
      //   sprintf(kernelName, "cdsScalarHaloExtract");
      //   cds->scalarHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
        
      //   sprintf(kernelName, "cdsScalarHaloScatter");
      //   cds->scalarHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);    
      // } 


      sprintf(kernelName, "cdsHaloGet");
      cds->haloGetKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
      
      sprintf(kernelName, "cdsHaloPut");
      cds->haloPutKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      if(cds->Nsubsteps){
        sprintf(kernelName, "cdsScalarHaloGet");
        cds->scalarHaloGetKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
        
        sprintf(kernelName, "cdsScalarHaloPut");
        cds->scalarHaloPutKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);    
      } 

      sprintf(fileName, DCDS "/okl/cdsAdvection%s.okl", suffix);

      // needed to be implemented
      sprintf(kernelName, "cdsAdvectionCubatureVolume%s", suffix);
      cds->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsAdvectionCubatureSurface%s", suffix);
      cds->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(kernelName, "cdsAdvectionVolume%s", suffix);
      cds->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsAdvectionSurface%s", suffix);
      cds->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      // ===========================================================================
      
      sprintf(fileName, DCDS "/okl/cdsHelmholtzRhs%s.okl", suffix);
      sprintf(kernelName, "cdsHelmholtzRhsEXTBDF%s", suffix);
      cds->helmholtzRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(fileName, DCDS "/okl/cdsHelmholtzBC%s.okl", suffix);
      sprintf(kernelName, "cdsHelmholtzIpdgBC%s", suffix);
      cds->helmholtzRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsHelmholtzBC%s", suffix);
      cds->helmholtzRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsHelmholtzAddBC%s", suffix);
      cds->helmholtzAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

    
      if(cds->Nsubsteps){
        // Note that resU and resV can be replaced with already introduced buffer
        cds->o_Ue     = mesh->device.malloc(cds->NVfields*Ntotal*sizeof(dfloat), cds->Ue);
        cds->o_Sd     = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->Sd);
        cds->o_resS   = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->resS);
        cds->o_rhsSd  = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->rhsSd);

        sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
        sprintf(kernelName, "scaledAddwOffset");
        cds->scaledAddKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(fileName, DCDS "/okl/cdsSubCycle%s.okl", suffix);
        sprintf(kernelName, "cdsSubCycleVolume%s", suffix);
        cds->subCycleVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "cdsSubCycleSurface%s", suffix);
        cds->subCycleSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "cdsSubCycleCubatureVolume%s", suffix);
        cds->subCycleCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "cdsSubCycleCubatureSurface%s", suffix);
        cds->subCycleCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(fileName, DCDS "/okl/cdsSubCycle.okl");
        sprintf(kernelName, "cdsSubCycleRKUpdate");
        cds->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "cdsSubCycleExt");
        cds->subCycleExtKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
   
    }
    MPI_Barrier(mesh->comm);
  }



  
#if 0
  printf(" Solver Parameters......\n");
  printf("alfa\t:\t %.8e \n", cds->alf);
  printf("invalfa\t:\t %.8e \n", cds->ialf);
  printf("k\t:\t %.8e \n", cds->k);
  printf("cp\t:\t %.8e \n", cds->cp);
  printf("rho\t:\t %.8e \n", cds->rho);
  printf("nu\t:\t %.8e \n", cds->nu);
  printf("Re\t:\t %.8e \n", cds->Re);
  printf("Pr\t:\t %.8e \n", cds->Pr);
  printf("dt\t:\t %.8e \n", cds->dt);
  printf("sdt\t:\t %.8e \n", cds->sdt);
  printf("invdt\t:\t %.8e \n", cds->idt);
  printf("ialf\t:\t %.8e \n", cds->ialf);
  printf("Nsubsteps\t:\t %02d \n", cds->Nsubsteps);
  printf("Nstages\t:\t %02d \n", cds->Nstages);
  printf("SNrk\t:\t %02d \n", cds->SNrk); 
  printf("ExplicitOrder\t:\t %02d \n", cds->ExplicitOrder);       
  
  cout<<kernelInfo;

#endif
  return cds;
}







