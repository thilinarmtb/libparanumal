#include "ins.h"

void insRunARK(ins_t *ins){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh_t *mesh = ins->mesh;
  
  //PID Control
  dfloat safe  = 0.8;   //safety factor

  //error control parameters
  dfloat beta    = 0.05; 
  dfloat factor1 = 0.333;
  dfloat factor2 = 3.0;
  dfloat exp1       = 0.25 - 0.75*beta;
  dfloat invfactor1 = 1.0/factor1;
  dfloat invfactor2 = 1.0/factor2;
  dfloat facold     = 1E-4;


  occa::initTimer(mesh->device);
  occaTimerTic(mesh->device,"INS");

  // Write Initial Data
  insReport(ins, 0.0, 0);

  int Nregect = 0;
  int allStep = 0;

  ins->tstep = 0;
  int done = 0;
  ins->time = ins->startTime;
  
  while (!done) {

    if (ins->dt<ins->dtMIN){
      printf("ERROR: Time step became too small at time step=%d\n", ins->tstep);
      exit (-1);
    }
    if (isnan(ins->dt)) {
      printf("ERROR: Solution became unstable at time step=%d\n", ins->tstep);
      exit (-1);
    }

    //check for final timestep
    if (ins->time+ins->dt > ins->finalTime){
      ins->dt = ins->finalTime-ins->time;
      done = 1;
    }

    insAdvection(ins, ins->time, ins->o_U, ins->o_NU);
    insDiffusion(ins, ins->time, ins->o_U, ins->o_LU);
    insGradient (ins, ins->time, ins->o_P, ins->o_GP);

    for(int stage=1;stage<=ins->Nstages;++stage){

      // intermediate stage time
      dfloat stageTime = ins->time + ins->rkC[stage]*ins->dt;

      insVelocityRhs  (ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
      insVelocitySolve(ins, stageTime, stage, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

      insPressureRhs  (ins, stageTime, stage);
      insPressureSolve(ins, stageTime, stage);      

      insPressureUpdate(ins, stageTime, stage, ins->o_rkP);
      insGradient(ins, stageTime, ins->o_rkP, ins->o_rkGP);

      insVelocityUpdate(ins, stageTime, stage, ins->o_rkGP, ins->o_rkU);
      
      //compute and save NU and LU
      insAdvection(ins, stageTime, ins->o_rkU, ins->o_rkNU);
      insDiffusion(ins, stageTime, ins->o_rkU, ins->o_rkLU); 


      ins->o_NU.copyFrom(ins->o_rkNU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_LU.copyFrom(ins->o_rkLU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), stage*ins->Ntotal*sizeof(dfloat), 0);
      ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_GP.copyFrom(ins->o_rkGP, ins->Ntotal*ins->NVfields*sizeof(dfloat), stage*ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
    } 
    
    if (ins->embeddedRKFlag==0) {//check if an embedded rk method is being used
      //accept the step and proceed
      ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
      ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), 0);
      ins->tstep++;
      ins->time += ins->dt;
    } else {

      // final stage
      dfloat stageTime = ins->time + ins->dt;
      int stage = ins->Nstages+1;

      //update and get error estimate
      dfloat err = insVelocityRkUpdate(ins, stageTime, ins->o_rkU);

      dfloat fac1 = pow(err,exp1);
      dfloat fac  = fac1/pow(facold,beta);

      fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
      dfloat dtnew = ins->dt/fac;

      if(err<1.0){

        insPressureRhs  (ins, stageTime, stage);
        insPressureSolve(ins, stageTime, stage);      

        insPressureUpdate(ins, stageTime, stage, ins->o_rkP);
        insGradient(ins, stageTime, ins->o_rkP, ins->o_rkGP);

        insVelocityUpdate(ins, stageTime, stage, ins->o_rkGP, ins->o_rkU);

        //accept the step and proceed
        ins->o_U.copyFrom(ins->o_rkU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 0);
        ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat), 0);
        ins->tstep++;
        ins->time += ins->dt;

        facold = mymax(err,1E-4);

      }else{
        Nregect++;
        dtnew = ins->dt/(mymax(invfactor1,fac1/safe));
        done =0;
      }

      ins->dt = dtnew;
      allStep++;
    }

    printf("\rTime = %.4e (%d). Average Dt = %.4e, Rejection rate = %.2g   \n", ins->time, ins->tstep, ins->time/(dfloat)ins->tstep, Nregect/(dfloat) ins->tstep); fflush(stdout);
    occaTimerTic(mesh->device,"Report");
    if(((ins->tstep)%(ins->outputStep))==0){
      if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterP);
      if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
      insReport(ins, ins->time, ins->tstep);
    }

    if (ins->dim==2 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d \n", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    if (ins->dim==3 && rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", ins->tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    occaTimerToc(mesh->device,"Report");


/*
    
  */
  }
  occaTimerToc(mesh->device,"INS");


  dfloat finalTime = ins->NtimeSteps*ins->dt;
  printf("\n");
  insReport(ins, finalTime, ins->NtimeSteps);
  
  if(rank==0) occa::printTimer();
}
