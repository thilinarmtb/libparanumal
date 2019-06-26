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

#include "ins.h"

void insRunEXTBDF(ins_t *ins){

  mesh_t *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  int NekSubCycle = 0;
  if(ins->options.compareArgs("SUBCYCLING TYPE", "NEK"))
    NekSubCycle = 1;
  
  int NstokesSteps = 0;
  dfloat oldDt = ins->dt;
  ins->dt *= 100;

  if (mesh->rank==0) printf("Number of Timesteps: %d\n", ins->NtimeSteps);
  
  for(int tstep=0;tstep<NstokesSteps;++tstep){
    if(tstep<1) 
      insExtBdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      insExtBdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3)
      insExtBdfCoefficents(ins,tstep+1);

    insGradient (ins, 0, ins->o_P, ins->o_GP);

    insVelocityRhs  (ins, 0, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, 0, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    insPressureRhs  (ins, 0, ins->Nstages);
    insPressureSolve(ins, 0, ins->Nstages); 

    insPressureUpdate(ins, 0, ins->Nstages, ins->o_rkP);
    insGradient(ins, 0, ins->o_rkP, ins->o_rkGP);

    //cycle history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_P.copyFrom(ins->o_P, ins->Ntotal*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*sizeof(dfloat));
    }

    //copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    //update velocity
    insVelocityUpdate(ins, 0, ins->Nstages, ins->o_rkGP, ins->o_rkU);

    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    if (mesh->rank==0) printf("\rSstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
  }

  if (mesh->rank==0) printf("\n");

  ins->dt = oldDt;
  // Write Initial Data
  if(ins->outputStep) insReport(ins, ins->startTime, 0);

  // for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
  for(int tstep=0;tstep<5;++tstep){

    // if(ins->restartedFromFile){
      // if(tstep=0 && ins->temporalOrder>=2) 
      //   insExtBdfCoefficents(ins,2);
      // else if(tstep=0 && ins->temporalOrder>=3) 
      //   insExtBdfCoefficents(ins,3);
    // }else{
      if(tstep<1) 
        insExtBdfCoefficents(ins,tstep+1);
      else if(tstep<2 && ins->temporalOrder>=2) 
        insExtBdfCoefficents(ins,tstep+1);
      else if(tstep<3 && ins->temporalOrder>=3) 
        insExtBdfCoefficents(ins,tstep+1);
    // }

    dfloat time = ins->startTime + tstep*ins->dt;

    dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

    if(ins->Nsubsteps) {
      if(!NekSubCycle)
	insSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
      else
	insNekSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
    } else {
      insAdvection(ins, time, ins->o_U, ins->o_NU);
    }
    
    insGradient (ins, time, ins->o_P, ins->o_GP);
    
    // Add explicit contrubitions like gravity, relaxation etc...
    insAddVelocityRhs(ins, time); 

    
    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    insPressureSolve(ins, time+ins->dt, ins->Nstages); 

    insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
    insGradient(ins, time+ins->dt, ins->o_rkP, ins->o_rkGP);

    //cycle history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_P.copyFrom(ins->o_P, ins->Ntotal*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*sizeof(dfloat));
    }

    //copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    //update velocity
    insVelocityUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkGP, ins->o_rkU);

    if(mesh->dim==3){
      if(ins->elementType == QUADRILATERALS ||
	 ins->elementType == TRIANGLES){
	ins->constrainKernel(mesh->Nelements,
			     offset,
			     mesh->o_x,
			     mesh->o_y,
			     mesh->o_z,
			     ins->o_rkU);
      }
    }
          
    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_FU.copyFrom(ins->o_FU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }



#if 0
    if(((tstep+1)%10)==0){
      char fname[BUFSIZ];
      sprintf(fname, "Iterations_Quad3D_%d", ins->SNrk);
      FILE *fp; 
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d %d %d %d\n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
      fclose(fp);
    }
#endif

    occaTimerTic(mesh->device,"Report");

    if(ins->outputStep){
      if(((tstep+1)%(ins->outputStep))==0){
        if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
        if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);

	insReport(ins, time+ins->dt, tstep+1);

        // Write a restart file
        if(ins->writeRestartFile){
          if(mesh->rank==0) printf("\nWriting Binary Restart File....");
            insRestartWrite(ins, ins->options, time+ins->dt);
          if(mesh->rank==0) printf("done\n");
        }

        // // Update Time-Step Size
        // if(ins->dtAdaptStep){
        //   if(((ins->tstep)%(ins->dtAdaptStep))==0){
        //     if(rank==0) printf("\n Adapting time Step Size to ");
        //       insComputeDt(ins, ins->time);
        //     if(rank==0) printf("%.4e\n", ins->dt);
        //      // Interpolate history for the new time step size
        //       insInterpolateHistory(ins, ins->dtold, ins->dt);

        //   }
        // } 
      }
    }

    if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
    occaTimerToc(mesh->device,"Report");
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");

  if(ins->outputStep) insReport(ins, finalTime,ins->NtimeSteps);
  
  if(mesh->rank==0) occa::printTimer();
}

