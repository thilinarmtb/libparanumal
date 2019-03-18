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

#include "stokes.h"

static void stokesJacobiPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv);
static void stokesSchurComplementBlockDiagPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv);

void stokesPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokesVecCopy(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    stokesJacobiPreconditioner(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "SCHURCOMPLEMENTBLOCKDIAG")) {
    stokesSchurComplementBlockDiagPreconditioner(stokes, v, Mv);
  } else {
    printf("ERROR:  Invalid value %s for [PRECONDITIONER] option.",
           stokes->options.getArgs("PRECONDITIONER").c_str());
  }

  return;
}

static void stokesJacobiPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv)
{
  stokes->dotMultiplyKernel(stokes->Ndof,
                            v.o_v,
                            stokes->precon->invDiagA.o_v,
                            Mv.o_v);
  return;
}

static void stokesSchurComplementBlockDiagPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv)
{
  stokesVec_t tmp, tmp2, tmp3;

  stokesVecAllocate(stokes, &tmp);
  stokesVecAllocate(stokes, &tmp2);
  stokesVecAllocate(stokes, &tmp3);

#if 1
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_x, Mv.o_x);
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_y, Mv.o_y);
  if (stokes->meshV->dim == 3)
    ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_z, Mv.o_z);

  stokes->dotMultiplyKernel(stokes->NtotalP, stokes->precon->invMM.o_p, v.o_p, Mv.o_p);

  stokes->dotMultiplyKernel(stokes->NtotalV, stokes->ogs->o_invDegree, Mv.o_x, Mv.o_x);
  stokes->dotMultiplyKernel(stokes->NtotalV, stokes->ogs->o_invDegree, Mv.o_y, Mv.o_y);
  if (stokes->meshV->dim == 3)
    stokes->dotMultiplyKernel(stokes->NtotalV, stokes->ogs->o_invDegree, Mv.o_z, Mv.o_z);
  stokes->dotMultiplyKernel(stokes->NtotalP, stokes->meshP->ogs->o_invDegree, Mv.o_p, Mv.o_p);

  ogsGatherScatter(Mv.o_x, ogsDfloat, ogsAdd, stokes->ogs);
  ogsGatherScatter(Mv.o_y, ogsDfloat, ogsAdd, stokes->ogs);
  if (stokes->meshV->dim == 3)
    ogsGatherScatter(Mv.o_z, ogsDfloat, ogsAdd, stokes->ogs);
  ogsGatherScatter(Mv.o_p, ogsDfloat, ogsAdd, stokes->meshP->ogs);

  if (stokes->Nmasked) {
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Mv.o_x);
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Mv.o_y);
    if (stokes->meshV->dim == 3)
      stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Mv.o_z);
  }
#else

  // zero pressure
  stokesVecZero(stokes, Mv);
  
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_x, Mv.o_x);
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_y, Mv.o_y);

  stokesOperator(stokes, Mv, tmp);

  stokesVecScaledAdd(stokes, 1.0, v, -1.0, tmp);

  ellipticPreconditioner(stokes->precon->elliptic, 0.0, tmp.o_x, tmp2.o_x);
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, tmp.o_y, tmp2.o_y);

  //  dfloat relax = 0.222;
  dfloat relax = 0;
  stokesVecScaledAdd(stokes, relax, tmp2, 1.0, Mv);

  stokes->dotMultiplyKernel(stokes->NtotalV, stokes->meshV->ogs->o_invDegree, Mv.o_x, Mv.o_x);
  ogsGatherScatter(Mv.o_x, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  stokes->dotMultiplyKernel(stokes->NtotalV, stokes->meshV->ogs->o_invDegree, Mv.o_y, Mv.o_y);
  ogsGatherScatter(Mv.o_y, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  
  stokes->dotMultiplyKernel(stokes->NtotalP, stokes->precon->invMM.o_p, v.o_p, Mv.o_p);
  stokes->dotMultiplyKernel(stokes->NtotalP, stokes->meshP->ogs->o_invDegree, Mv.o_p, Mv.o_p);
  ogsGatherScatter(Mv.o_p, ogsDfloat, ogsAdd, stokes->meshP->ogs);
#endif

  stokesVecFree(stokes, &tmp);
  stokesVecFree(stokes, &tmp2);
  stokesVecFree(stokes, &tmp3);

  return;
}
