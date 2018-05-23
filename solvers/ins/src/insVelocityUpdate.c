#include "ins.h"

void insVelocityUpdate(ins_t *ins, dfloat time, int stage, 
                        occa::memory o_rkGP,
                        occa::memory o_rkU){

  mesh_t *mesh = ins->mesh;

  occa::memory o_prkA = ins->o_prkA;
  occa::memory o_prkAX = ins->o_prkAX;
  int finalStageSwitch = 0;
  if (stage==ins->Nstages+1) { //final stage
    stage = ins->Nstages;
    finalStageSwitch = 1;
    o_prkA = ins->o_prkB;
    o_prkAX = ins->o_prkBX;
  }

  // U^s = Uhat - dt * GP^s + dt*\sum^s-1 pa_si GP^i
  occaTimerTic(mesh->device,"VelocityUpdate");
  ins->velocityUpdateKernel(mesh->Nelements,
                              stage,
                              ins->ARKswitch,
                              finalStageSwitch,
                              ins->dt,
                              ins->fieldOffset,
                              o_prkA,
                              o_prkAX,
                              o_rkGP,
                              ins->o_GP,
                              o_rkU);
  occaTimerToc(mesh->device,"VelocityUpdate");
}


void insVelocityRkUpdate(ins_t *ins, dfloat time, occa::memory o_rkU){

  mesh_t *mesh = ins->mesh;

  // U^s = Uhat - dt * GP^s + dt*\sum^s-1 pa_si GP^i
  occaTimerTic(mesh->device,"VelocityUpdate");
  ins->velocityRkUpdateKernel(mesh->Nelements,
                              ins->Nstages+1,
                              ins->dt,
                              ins->fieldOffset,
                              ins->o_erkB,
                              ins->o_irkB,
                              ins->o_prkB,
                              ins->o_prkBX,
                              ins->o_U,
                              ins->o_NU,
                              ins->o_LU,
                              ins->o_GP,
                              o_rkU);
  occaTimerToc(mesh->device,"VelocityUpdate");
}

