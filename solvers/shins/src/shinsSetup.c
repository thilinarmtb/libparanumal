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

#include "shins.h"
#include "omp.h"
#include <unistd.h>

shins_t *shinsSetup(mesh_t *mesh, setupAide options){

  shins_t *shins = (shins_t*) calloc(1, sizeof(shins_t));
  shins->mesh = mesh;
  shins->options = options;

  options.getArgs("MESH DIMENSION", shins->dim);
  options.getArgs("ELEMENT TYPE", shins->elementType);

  shins->NVfields = (shins->dim==3) ? 3:2;  // Total Number of Velocity Fields
  shins->NTfields = (shins->dim==3) ? 4:3;  // Total Velocity + Pressure

  mesh->Nfields = 1; 

  shins->g0 =  1.0;

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    shins->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    shins->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    shins->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));

    shins->extC = (dfloat*) calloc(3, sizeof(dfloat));
  }else{ printf("Error: shins uses only EXBDF currently"); exit(-1); }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    shins->Nstages = 1;
    shins->temporalOrder = 1;
    shins->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    shins->Nstages = 2;
    shins->temporalOrder = 2;
    shins->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    shins->Nstages = 3;
    shins->temporalOrder = 3;
    shins->g0 = 11.f/6.f;
  }

  shins->readRestartFile = 0; 
  options.getArgs("RESTART FROM FILE", shins->readRestartFile);
  
  shins->writeRestartFile = 0; 
  options.getArgs("WRITE RESTART FILE", shins->writeRestartFile);

   // Set radial expansion
  shins->Nmodes = 0; 
  options.getArgs("RADIAL EXPANSION MODES", shins->Nmodes);

  shins->R = 1.0;
  options.getArgs("OUTER RADIUS", shins->R);


  shins->Nlocal = mesh->Np*mesh->Nelements;
  shins->Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  shins->fieldOffset = shins->Nmodes*shins->Ntotal;

  shins->Nblock      = (shins->Nlocal+blockSize-1)/blockSize; // Check this Nlocal size!!!!!AK.

  // compute samples of q at interpolation nodes
  shins->U     = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Nstages*shins->Ntotal,sizeof(dfloat));
  shins->P     = (dfloat*) calloc(                shins->Nmodes*shins->Nstages*shins->Ntotal,sizeof(dfloat));

  //rhs storage
  shins->rhsU  = (dfloat*) calloc(shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rhsV  = (dfloat*) calloc(shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rhsW  = (dfloat*) calloc(shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rhsP  = (dfloat*) calloc(shins->Nmodes*shins->Ntotal,sizeof(dfloat));

  //additional field storage
  shins->NU   = (dfloat*) calloc(shins->NVfields*(shins->Nstages+1)*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->LU   = (dfloat*) calloc(shins->NVfields*(shins->Nstages+1)*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->GP   = (dfloat*) calloc(shins->NVfields*(shins->Nstages+1)*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  
  shins->GU   = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal*4,sizeof(dfloat));
  
  shins->rkU  = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rkP  = (dfloat*) calloc(                shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->PI   = (dfloat*) calloc(                shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  
  shins->rkNU = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rkLU = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->rkGP = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));

  //plotting fields
  shins->Vort = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  shins->Div  = (dfloat*) calloc(                shins->Nmodes*shins->Nlocal,sizeof(dfloat));
  
  
  shins->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",shins->Nsubsteps);

  if(shins->Nsubsteps){
    shins->Ud    = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
    shins->Ue    = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
    shins->resU  = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
    shins->rhsUd = (dfloat*) calloc(shins->NVfields*shins->Nmodes*shins->Ntotal,sizeof(dfloat));
  }

  dfloat rho  = 1.0 ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", shins->ubar);
  options.getArgs("VBAR", shins->vbar);
  if (shins->dim==3)
    options.getArgs("WBAR", shins->wbar);
  options.getArgs("PBAR", shins->pbar);
  options.getArgs("VISCOSITY", shins->nu);

  //Reynolds number
  shins->Re = shins->ubar/shins->nu;

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(shins->elementType==QUADRILATERALS)
    meshOccaSetupQuad3D(mesh, options, kernelInfo); 
  else if(shins->elementType==TRIANGLES){
    printf("Error OccaSetup: Only Quad3D is supported currently\n");exit(-1);
    // meshOccaSetupTri3D(mesh, options, kernelInfo);
  }
  else{ printf("Error OccaSetup: Only Quad3D and Tri3D are suppurted\n");exit(-1);}  

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_pbar"]= shins->pbar;
  kernelInfo["defines/" "p_ubar"]= shins->ubar;
  kernelInfo["defines/" "p_vbar"]= shins->vbar;
  kernelInfo["defines/" "p_wbar"]= shins->wbar;
  kernelInfo["defines/" "p_nu"]  = shins->nu;

  kernelInfo["defines/" "p_NTfields"]= shins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= shins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  shins->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  shins->Nsubsteps;

  kernelInfo["defines/" "p_fainv"] = (dfloat) 0.0;
  kernelInfo["defines/" "p_invRadiusSq"] = (dfloat) 1./(mesh->sphereRadius*mesh->sphereRadius);

  // AK: No need to hold it but may be needed later // add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  shins->o_U = mesh->device.malloc(shins->NVfields*shins->Nstages*shins->Nmodes*shins->Ntotal*sizeof(dfloat), shins->U);
  shins->o_P = mesh->device.malloc(                shins->Nstages*shins->Nmodes*shins->Ntotal*sizeof(dfloat), shins->P);

#if 0
  if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);
#endif
  
// Assume Zero initial conditions currently, !!! modify later
#if 0
if(options.compareArgs("INITIAL CONDITION", "BROWN-MINION") &&
  (ins->elementType == QUADRILATERALS && ins->dim==3)){
  printf("Setting up initial condition for BROWN-MINION test case...");
  insBrownMinionQuad3D(ins);
  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);
  printf("done\n");

  char fname[BUFSIZ];
  string outName;
  ins->options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, ins->frame++);
  insPlotVTU(ins, fname);
}else{

  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {
      if (ins->dim==2) 
        ins->setFlowFieldKernel =  mesh->device.buildKernel(DINS "/okl/insSetFlowField2D.okl", "insSetFlowField2D", kernelInfo);  
      else
        ins->setFlowFieldKernel =  mesh->device.buildKernel(DINS "/okl/insSetFlowField3D.okl", "insSetFlowField3D", kernelInfo);  
    }
    MPI_Barrier(mesh->comm);
  }

  ins->startTime =0.0;
  options.getArgs("START TIME", ins->startTime);
  ins->setFlowFieldKernel(mesh->Nelements,
                          ins->startTime,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          ins->fieldOffset,
                          ins->o_U,
                          ins->o_P);
  ins->o_U.copyTo(ins->U);

}

#endif

  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }

    for(int n=0;n<mesh->Np;++n){
      for(int m=0; m<shins->Nmodes; m++){
        const dlong id = n + mesh->Np*e + m*mesh->Np*mesh->Nelements;
        dfloat t = 0;
        dfloat uxn = shins->U[id+0*shins->fieldOffset]; // this time -> fieldOffset = Nmodes*Ntotal(e*Np)
        dfloat uyn = shins->U[id+1*shins->fieldOffset]; // this time -> fieldOffset = Nmodes*Ntotal(e*Np)
        dfloat uzn = shins->U[id+2*shins->fieldOffset]; // this time -> fieldOffset = Nmodes*Ntotal(e*Np)
        //Squared maximum velocity
        dfloat numax = uxn*uxn + uyn*uyn + uzn*uzn;
        umax = mymax(umax, numax); 
      }
    }
  }


  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity

  options.getArgs("CFL", shins->cfl);
  dfloat dt     = shins->cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(shins->dti), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  // save initial time-step estimate 
  shins->dt = shins->dti; 

  options.getArgs("FINAL TIME", shins->finalTime);
  options.getArgs("START TIME", shins->startTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    shins->NtimeSteps = (shins->finalTime-shins->startTime)/shins->dt;

    if(shins->Nsubsteps){
      shins->dt         = shins->Nsubsteps*shins->dt;
      shins->NtimeSteps = (shins->finalTime-shins->startTime)/shins->dt;
      shins->dt         = (shins->finalTime-shins->startTime)/shins->NtimeSteps;
      shins->sdt        = shins->dt/shins->Nsubsteps;
    } else{
      shins->NtimeSteps = (shins->finalTime-shins->startTime)/shins->dt;
      shins->dt         = (shins->finalTime-shins->startTime)/shins->NtimeSteps;
    }
  }

  shins->dtMIN = 1E-2*shins->dt; //minumum allowed timestep

  if (mesh->rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  shins->cfl);
    printf("dt = %g\n",   dt);
  }

  if (shins->Nsubsteps && mesh->rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", shins->dt, shins->sdt, shins->dt/shins->sdt);
  
  // Hold some inverses for kernels
  shins->inu    = 1.0/shins->nu; 
  shins->idt    = 1.0/shins->dt;
  shins->lambda = shins->g0 / (shins->dt * shins->nu);


  //make option objects for elliptc solvers
  shins->vOptions = options;
  shins->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  shins->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  shins->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  shins->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  shins->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  shins->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  shins->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  shins->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  shins->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));

  shins->pOptions = options;
  shins->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  shins->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  shins->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  shins->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  shins->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  shins->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  shins->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  shins->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  shins->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));

  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  // AK: No need to BC contions currently 
  int uBCType[7] = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int vBCType[7] = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int wBCType[7] = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[7] = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances 
  shins->presTOL = 1E-8;
  shins->velTOL  = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  shins->uSolver = (asbf_t*) calloc(1, sizeof(asbf_t));
  shins->uSolver->mesh = mesh;
  shins->uSolver->options = shins->vOptions;
  shins->uSolver->dim = shins->dim;
  shins->uSolver->elementType = shins->elementType;
  shins->uSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(shins->uSolver->BCType,uBCType,7*sizeof(int));
  asbfSolveSetup(shins->uSolver, shins->lambda, kernelInfoV); 

#if 0
  shins->vSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  shins->vSolver->mesh = mesh;
  shins->vSolver->options = shins->vOptions;
  shins->vSolver->dim = shins->dim;
  shins->vSolver->elementType = shins->elementType;
  shins->vSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(shins->vSolver->BCType,vBCType,7*sizeof(int));
  ellipticSolveSetup(shins->vSolver, shins->lambda, kernelInfoV); //!!!!!

  if (shins->dim==3) {
    shins->wSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
    shins->wSolver->mesh = mesh;
    shins->wSolver->options = shins->vOptions;
    shins->wSolver->dim = shins->dim;
    shins->wSolver->elementType = shins->elementType;
    shins->wSolver->BCType = (int*) calloc(7,sizeof(int));
    memcpy(shins->wSolver->BCType,wBCType,7*sizeof(int));
    ellipticSolveSetup(shins->wSolver, shins->lambda, kernelInfoV);  //!!!!! 
  }
  
  if (mesh->rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  shins->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  shins->pSolver->mesh = mesh;
  shins->pSolver->options = shins->pOptions;
  shins->pSolver->dim = shins->dim;
  shins->pSolver->elementType = shins->elementType;
  shins->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(shins->pSolver->BCType,pBCType,7*sizeof(int));
  ellipticSolveSetup(shins->pSolver, 0.0, kernelInfoP); //!!!!


  //make node-wise boundary flags
  shins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  shins->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) shins->VmapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          shins->VmapB[fid+e*mesh->Np] = mymin(bc,shins->VmapB[fid+e*mesh->Np]);
          shins->PmapB[fid+e*mesh->Np] = mymax(bc,shins->PmapB[fid+e*mesh->Np]);
        }
      }
    }
  }

  ogsGatherScatter(shins->VmapB, ogsInt, ogsMin, mesh->ogs);
  ogsGatherScatter(shins->PmapB, ogsInt, ogsMax, mesh->ogs);

  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (shins->VmapB[n] == 1E9) {
      shins->VmapB[n] = 0.;
    }
  }
  shins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), shins->VmapB);
  shins->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), shins->PmapB);


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

  
  // IsoSurface related
  if(ins->dim==3 && ins->elementType != QUADRILATERALS ){
    kernelInfo["defines/" "p_isoNfields"]= ins->isoNfields;
    // Define Isosurface Area Tolerance
    kernelInfo["defines/" "p_triAreaTol"]= (dfloat) 1.0E-16;

    kernelInfo["defines/" "p_dim"]= ins->dim;
    kernelInfo["defines/" "p_plotNp"]= mesh->plotNp;
    kernelInfo["defines/" "p_plotNelements"]= mesh->plotNelements;
    
    int plotNthreads = mymax(mesh->Np, mymax(mesh->plotNp, mesh->plotNelements));
    kernelInfo["defines/" "p_plotNthreads"]= plotNthreads;
 } 





  // if (mesh->rank==0) {
  //   printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  //   printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

  //   printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  // }
  
  if (options.compareArgs("TIME INTEGRATOR", "ARK")) {
    ins->o_rkC  = mesh->device.malloc(         ins->Nrk*sizeof(dfloat),ins->rkC );
    ins->o_erkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->erkA);
    ins->o_irkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->irkA);
    ins->o_prkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->prkA);
    ins->o_prkB = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->prkB);
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    ins->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_prkA = ins->o_extbdfC;
    ins->o_prkB = ins->o_extbdfC;
  }

  // MEMORY ALLOCATION
  ins->o_rhsU  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsP);

  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_LU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->LU);
  ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  
  ins->o_GU    = mesh->device.malloc(ins->NVfields*Ntotal*4*sizeof(dfloat), ins->GU);
  
  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);
  
  ins->o_rkNU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkNU);
  ins->o_rkLU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkLU);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);

  //storage for helmholtz solves
  ins->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  //plotting fields
  ins->o_Vort = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Vort);
  ins->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), ins->Div);

  if(ins->elementType==HEXAHEDRA) // !!!! check that
    ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);
  else 
    ins->o_cU = ins->o_U;

  if(mesh->totalHaloPairs){//halo setup
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes = ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);
    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);

#if 0
    occa::memory o_vsendBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    occa::memory o_vrecvBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    occa::memory o_psendBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    occa::memory o_precvBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(vGatherBytes, NULL);
    
    ins->vSendBuffer = (dfloat*) o_vsendBuffer.getMappedPointer();
    ins->vRecvBuffer = (dfloat*) o_vrecvBuffer.getMappedPointer();
    ins->pSendBuffer = (dfloat*) o_psendBuffer.getMappedPointer();
    ins->pRecvBuffer = (dfloat*) o_precvBuffer.getMappedPointer();
    ins->velocityHaloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
#endif
    occa::memory o_vSendBuffer, o_vRecvBuffer, o_pSendBuffer, o_pRecvBuffer, o_gatherTmpPinned;

    ins->vSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vSendBuffer);
    ins->vRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vRecvBuffer);

    ins->pSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pSendBuffer);
    ins->pRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pRecvBuffer);

    ins->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, ins->o_gatherTmpPinned);
    
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  ins->velocityHaloGatherTmp);
  }

  // set kernel name suffix
  char *suffix;
  
  if(ins->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(ins->elementType==QUADRILATERALS){
    if(ins->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D"); 
  }
  if(ins->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(ins->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {

      sprintf(fileName, DINS "/okl/insHaloExchange.okl");
      sprintf(kernelName, "insVelocityHaloExtract");
      ins->velocityHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insVelocityHaloScatter");
      ins->velocityHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloExtract");
      ins->pressureHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloScatter");
      ins->pressureHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // --
      if(ins->dim==3 && ins->elementType==QUADRILATERALS){
	sprintf(fileName, DINS "/okl/insConstrainQuad3D.okl");
	sprintf(kernelName, "insConstrainQuad3D");
	ins->constrainKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
      
      // ===========================================================================

      sprintf(fileName, DINS "/okl/insAdvection%s.okl", suffix);

      // needed to be implemented
      sprintf(kernelName, "insAdvectionCubatureVolume%s", suffix);
      ins->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionCubatureSurface%s", suffix);
      ins->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionVolume%s", suffix);
      ins->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionSurface%s", suffix);
      ins->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insDiffusion%s.okl", suffix);
      // sprintf(kernelName, "insDiffusion%s", suffix);
      // ins->diffusionKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insDiffusionIpdg%s.okl", suffix);
      // sprintf(kernelName, "insDiffusionIpdg%s", suffix);
      // ins->diffusionIpdgKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insVelocityGradient%s.okl", suffix);
      // sprintf(kernelName, "insVelocityGradient%s", suffix);
      // ins->velocityGradientKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // // ===========================================================================

      sprintf(fileName, DINS "/okl/insGradient%s.okl", suffix);
      sprintf(kernelName, "insGradientVolume%s", suffix);
      ins->gradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insGradientSurface%s", suffix);
      ins->gradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insDivergence%s.okl", suffix);
      sprintf(kernelName, "insDivergenceVolume%s", suffix);
      ins->divergenceVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insDivergenceSurface%s", suffix);
      ins->divergenceSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insVelocityRhs%s.okl", suffix);
      if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
        sprintf(kernelName, "insVelocityRhsARK%s", suffix); 
      else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) 
        sprintf(kernelName, "insVelocityRhsEXTBDF%s", suffix);
      ins->velocityRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);



      if(!(ins->dim==3 && ins->elementType==QUADRILATERALS) ){
        sprintf(fileName, DINS "/okl/insVelocityBC%s.okl", suffix);
        sprintf(kernelName, "insVelocityIpdgBC%s", suffix);
        ins->velocityRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insVelocityBC%s", suffix);
        ins->velocityRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insVelocityAddBC%s", suffix);
        ins->velocityAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // ===========================================================================
      
      // Dont forget to modify!!!!!!!!
      sprintf(fileName, DINS "/okl/insPressureRhs%s.okl", suffix);
      sprintf(kernelName, "insPressureRhs%s", suffix);
      ins->pressureRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

       if(!(ins->dim==3 && ins->elementType==QUADRILATERALS) ){
        sprintf(fileName, DINS "/okl/insPressureBC%s.okl", suffix);
        sprintf(kernelName, "insPressureIpdgBC%s", suffix);
        ins->pressureRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insPressureBC%s", suffix);
        ins->pressureRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insPressureAddBC%s", suffix);
        ins->pressureAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insPressureUpdate.okl");
      sprintf(kernelName, "insPressureUpdate");
      ins->pressureUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insVelocityUpdate.okl");
      sprintf(kernelName, "insVelocityUpdate");
      ins->velocityUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);      

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insVorticity%s.okl", suffix);
      sprintf(kernelName, "insVorticity%s", suffix);
      ins->vorticityKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
    
      // ===========================================================================
      if(ins->dim==3 && ins->options.compareArgs("OUTPUT TYPE","ISO")){
        sprintf(fileName, DINS "/okl/insIsoSurface3D.okl");
        sprintf(kernelName, "insIsoSurface3D");

        ins->isoSurfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
      }
      

      // Not implemented for Quad 3D yet !!!!!!!!!!
      if(ins->Nsubsteps){
        // Note that resU and resV can be replaced with already introduced buffer
        ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
        ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
        ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
        ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

        if(ins->elementType==HEXAHEDRA)
          ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
        else 
          ins->o_cUd = ins->o_Ud;

        sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
        sprintf(kernelName, "scaledAddwOffset");
        ins->scaledAddKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(fileName, DINS "/okl/insSubCycle%s.okl", suffix);
        sprintf(kernelName, "insSubCycleVolume%s", suffix);
        ins->subCycleVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleSurface%s", suffix);
        ins->subCycleSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleCubatureVolume%s", suffix);
        ins->subCycleCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleCubatureSurface%s", suffix);
        ins->subCycleCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(fileName, DINS "/okl/insSubCycle.okl");
        sprintf(kernelName, "insSubCycleRKUpdate");
        ins->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        sprintf(kernelName, "insSubCycleExt");
        ins->subCycleExtKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
    }
    MPI_Barrier(mesh->comm);
  }
#endif
  return shins;
}







