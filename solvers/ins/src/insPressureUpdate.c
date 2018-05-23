#include "ins.h"

void insPressureUpdate(ins_t *ins, dfloat time, int stage, occa::memory o_rkP){

  mesh_t *mesh = ins->mesh;
  
  occa::memory o_prkAX = ins->o_prkAX;
  int finalStageSwitch = 0;
  if (stage==ins->Nstages+1) { //final stage
    stage = ins->Nstages;
    finalStageSwitch = 1;
    o_prkAX = ins->o_prkBX;
  }

  // P^s = PI + \sum^s-1 prk_si P^i 
  occaTimerTic(mesh->device,"PressureUpdate");
  ins->pressureUpdateKernel(mesh->Nelements,
                            stage,
                            finalStageSwitch,
		                        o_prkAX,
                            ins->fieldOffset,
                            ins->o_PI,
                            ins->o_P,
                            o_rkP);
  occaTimerToc(mesh->device,"PressureUpdate");
}
