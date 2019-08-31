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

void insRunTOMBO(ins_t *ins){

  mesh_t *mesh = ins->mesh;
  
  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  int StrongSubCycle = 0;
  if(ins->options.compareArgs("ADVECTION TYPE", "CONVECTIVE"))
    StrongSubCycle = 1;

  ins->frame=1;  
  // Write Initial Data
  if(ins->outputStep) insReport(ins, ins->startTime, 0);

  dfloat cfl = insComputeCfl(ins, ins->startTime, 0); printf("CFL = %.4e \n", cfl);
  
  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){
    // for(int tstep=0;tstep<10;++tstep){
    if(tstep<1) 
      insExtBdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      insExtBdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      insExtBdfCoefficents(ins,tstep+1);

    dfloat time = ins->startTime + tstep*ins->dt;

    dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
    if(ins->Nsubsteps) {
      if(!StrongSubCycle){
        insSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
      }else{
        insStrongSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
      }
    } else {
      insAdvection(ins, time, ins->o_U, ins->o_NU);
    }
    
  
    insCurlCurl(ins, time, ins->o_U, ins->o_NC); 

    // Add explicit contrubitions like gravity, relaxation etc...
    insAddVelocityRhs(ins, time); 

    insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    insPressureSolve(ins, time+ins->dt, ins->Nstages); 
    insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
   

    // copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 


    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    //cycle velocity history after velocityUpdate
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			(s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			(s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

   
    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 


    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_NC.copyFrom(ins->o_NC, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_FU.copyFrom(ins->o_FU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    

    occaTimerTic(mesh->device,"Report");

    if(ins->outputStep){
      if(((tstep+1)%(ins->outputStep))==0){
        if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
        if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);

	insReport(ins, time+ins->dt, tstep+1);
        
        dfloat cfl = insComputeCfl(ins, time+ins->dt, tstep+1); printf("CFL = %.4e \n", cfl);
    

        // Write a restart file
        if(ins->writeRestartFile){
          if(mesh->rank==0) printf("\nWriting Binary Restart File....");
	  insRestartWrite(ins, ins->options, time+ins->dt);
          if(mesh->rank==0) printf("done\n");
        }
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


