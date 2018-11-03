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

// solve lambda*U + A*U = rhsU
void insVelocitySolve(ins_t *ins, dfloat time, int stage, int velId, occa::memory o_Uhat){
  
  mesh_t *mesh = ins->mesh;
  elliptic_t *solver;
  if(velId==0) 
    solver = ins->uSolver; 
  else if(velId==1)
    solver = ins->vSolver; 
  else if(velId==2)
    solver = ins->wSolver;

  int quad3D = (ins->dim==3 && ins->elementType==QUADRILATERALS) ? 1 : 0;  
  
  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")){

    if(!quad3D) 
      ins->velocityRhsBCKernel(mesh->Nelements,
                                mesh->o_ggeo,
                                mesh->o_sgeo,
                                mesh->o_Dmatrices,
                                mesh->o_Smatrices,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                velId,
    	                          mesh->o_EToB,
                                mesh->o_sMT,
                                ins->lambda,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                ins->o_VmapB,
                                ins->o_rhsTmp);

    // gather-scatter
    ogsGatherScatter(ins->rhsTmp, ogsDfloat, ogsAdd, mesh->ogs);
    
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsTmp);
    
  } else if (ins->vOptions.compareArgs("DISCRETIZATION","IPDG") && !quad3D) {

    occaTimerTic(mesh->device,"velocityRhsIpdg");    
    ins->velocityRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  velId,
                                  solver->tau,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsTmp);
    occaTimerToc(mesh->device,"velocityRhsIpdg");   
  }

  //copy current velocity fields as initial guess? (could use Uhat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,velId*ins->fieldOffset*sizeof(dfloat));

  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_UH);
  }
  
  int Niter = 0; 
  occaTimerTic(mesh->device,"Ux-Solve");
  Niter = ellipticSolve(solver, ins->lambda, ins->velTOL, ins->o_rhsTmp, ins->o_UH);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  if(velId==0) ins->NiterU = Niter;
  if(velId==1) ins->NiterV = Niter;
  if(velId==2) ins->NiterW = Niter;

  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    ins->velocityAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            velId,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_VmapB,
                            ins->o_UH);
  }

  //copy into intermediate stage storage
  ins->o_UH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
}
