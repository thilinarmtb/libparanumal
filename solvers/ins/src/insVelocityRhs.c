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

void insVelocityRhs(ins_t *ins, dfloat time, int stage, occa::memory o_rhsU, occa::memory o_rhsV, occa::memory o_rhsW){
  
  mesh_t *mesh = ins->mesh; 

  if (ins->options.compareArgs("TIME INTEGRATOR", "ARK")) { 
    // rhsU^s = MM*(U^n - \sum^s-1 ea_si N(U^i) + \sum^s-1 ia_si LU^i - \sum^s-1 pa_si GP^i)/ia_ss nu dt
    ins->velocityRhsKernel(mesh->Nelements,
                           stage,
                           mesh->o_vgeo,
                           mesh->o_MM,
                           ins->idt,
                           ins->inu,
                           ins->o_erkA,
                           ins->o_irkA,
                           ins->o_prkA,
                           ins->o_prkB,
                           ins->fieldOffset,
                           ins->o_U,
                           ins->o_NU,
                           ins->o_LU,
                           ins->o_GP,
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);
  } else if (ins->options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {

    // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
    ins->velocityRhsKernel(mesh->Nelements,
                           mesh->o_vgeo,
                           mesh->o_MM,
                           ins->idt,
                           ins->inu,
                           ins->o_extbdfA,
                           ins->o_extbdfB,
                           ins->o_extbdfC,
                           ins->fieldOffset,
                           ins->o_U,
                           ins->o_NU,
                           ins->o_GP,
                           ins->o_FU,
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);
  }else if (ins->options.compareArgs("TIME INTEGRATOR", "TOMBO")) {

    // Get ready to velocity solve, compute Gradient of Pressure
    insGradient (ins, time, ins->o_P, ins->o_rkGP); 

    // rhsU^s = MM*1/nu*[ -(grad P) + sum_i ( (a_i) U^(n-q)/dt - b_i NU ^(n-q)) 
    ins->velocityRhsKernel(mesh->Nelements,
                           mesh->o_vgeo,
                           mesh->o_MM,
                           ins->idt,
                           ins->inu,
                           ins->g0,
                           ins->o_extbdfA,
                           ins->o_extbdfB,
                           ins->fieldOffset,
                           ins->o_U,
                           ins->o_NU,
                           ins->o_FU,
                           ins->o_rkGP,
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);

    // for Velocity Increament Solve
    if (ins->options.compareArgs("TIME INTEGRATOR", "DIFF")){
    
      ins->o_rkU.copyFrom(ins->o_U,ins->NVfields*ins->Ntotal*sizeof(dfloat),0,0);

      // Update Un+1 at dirichlet boundaries only
      ins->velocityAddBCTOMBOKernel(mesh->Nelements,
            time,
            ins->fieldOffset,
            mesh->o_sgeo,
            mesh->o_x,
            mesh->o_y,
            mesh->o_z,
            mesh->o_vmapM,
            ins->o_VmapB,   
            ins->o_rkU);

      // add [lap(Un) + lambda* Un] to Rhs to convert to velocity Difference
      //      const dfloat lambda = -ins->g0*ins->idt*ins->inu;  
      const dfloat lambda = -ins->lambda;  
      // simple AX kernel, modify later!!!!!!!AK 
      dlong voffset = 0*ins->fieldOffset; 
      // First for U
      ins->velocityAxKernel(mesh->Nelements,
          mesh->o_ggeo, 
          mesh->o_Dmatrices, 
          mesh->o_Smatrices, 
          mesh->o_MM, 
          lambda, 
          voffset, 
          ins->o_rkU, 
          o_rhsU); 

      voffset = 1*ins->fieldOffset; 
      // First for U
      ins->velocityAxKernel(mesh->Nelements,
          mesh->o_ggeo, 
          mesh->o_Dmatrices, 
          mesh->o_Smatrices, 
          mesh->o_MM, 
          lambda, 
          voffset, 
          ins->o_rkU, 
          o_rhsV); 


      if(ins->dim==3){
        voffset = 2*ins->fieldOffset; 
        ins->velocityAxKernel(mesh->Nelements,
                  mesh->o_ggeo, 
                  mesh->o_Dmatrices, 
                  mesh->o_Smatrices, 
                  mesh->o_MM, 
                  lambda, 
                  voffset, 
                  ins->o_rkU, 
                  o_rhsW); 
      }      
    }

  }
}