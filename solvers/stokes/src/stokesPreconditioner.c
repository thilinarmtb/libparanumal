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
static void stokesSCBlockPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv);

void stokesPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokesVecCopy(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    stokesJacobiPreconditioner(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "SCBLOCK")) {
    stokesSCBlockPreconditioner(stokes, v, Mv);
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

static void stokesSCBlockPreconditioner(stokes_t *stokes, stokesVec_t v, stokesVec_t Mv)
{

#if 0
  ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_x, Mv.o_x);
  ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_y, Mv.o_y);
  if (stokes->mesh->dim == 3)
    ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_z, Mv.o_z);
#else
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_x, Mv.o_x);
  ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_y, Mv.o_y);
  if (stokes->mesh->dim == 3)
    ellipticPreconditioner(stokes->precon->elliptic, 0.0, v.o_z, Mv.o_z);
#endif
 
  Mv.o_p.copyFrom(v.o_p, stokes->Ntotal*sizeof(dfloat));

  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
                                  stokes->o_vP,
                                  stokes->o_uP,
                                  Mv.o_p,
				  Mv.o_p);

  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
				  stokes->o_uP,
                                  stokes->o_vP,
                                  Mv.o_p,
				  Mv.o_p);
  
  return;

#if 0
  stokesVecZero(stokes, v);
  stokesVecZero(stokes, Mv);

  stokesVec_t e;

  stokesVecAllocate(stokes, &e);
  for (int i = 0; i < stokes->Ntotal; i++) {
    e.x[i] = 1.0;
    e.y[i] = 1.0;
    e.p[i] = 0.0;
  }

  stokesVecCopyHostToDevice(e);
  stokesVecGatherScatter(stokes, e);

  if (stokes->Nmasked) {
    stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, e.o_x);
    stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, e.o_y);
    if (stokes->mesh->dim == 3)
      stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, e.o_z);
  }

  stokesVecCopyDeviceToHost(e);
  printf("e:\n");
  stokesVecPrint(stokes, e);

  stokesOperator(stokes, e, v);

  stokesVecCopyDeviceToHost(v);
  printf("Ae:\n");
  stokesVecPrint(stokes, v);

  ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_x, Mv.o_x);
  ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_y, Mv.o_y);
  if (stokes->mesh->dim == 3)
    ellipticSolve(stokes->precon->elliptic, 0.0, 1.0e-8, v.o_z, Mv.o_z);

  Mv.o_p.copyFrom(v.o_p, stokes->Ntotal*sizeof(dfloat));

  stokesVecCopyDeviceToHost(Mv);
  printf("MAe:\n");
  stokesVecPrint(stokes, Mv);

  dfloat err = 0.0;
  for (int i = 0; i < stokes->Ntotal; i++) {
    err += pow(e.x[i] - Mv.x[i], 2.0);
    err += pow(e.y[i] - Mv.y[i], 2.0);
  }
  err = sqrt(err);

  printf("err = % .15e\n", err);

  exit(1);
#endif
}
