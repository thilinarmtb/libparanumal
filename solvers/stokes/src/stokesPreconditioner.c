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

static void stokesJacobiPreconditioner(stokes_t *stokes, occa::memory &v, occa::memory &Mv);
static void stokesSchurComplementBlockDiagPreconditioner(stokes_t *stokes, dfloat lambda, occa::memory &v, occa::memory &Mv);

void stokesPreconditioner(stokes_t *stokes, dfloat lambda,
			  occa::memory &v, occa::memory &Mv)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokesVecCopy(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    stokesJacobiPreconditioner(stokes, v, Mv);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "SCHURCOMPLEMENTBLOCKDIAG")) {
    stokesSchurComplementBlockDiagPreconditioner(stokes, lambda, v, Mv);
  } else {
    printf("ERROR:  Invalid value %s for [PRECONDITIONER] option.",
           stokes->options.getArgs("PRECONDITIONER").c_str());
  }

  return;
}

static void stokesJacobiPreconditioner(stokes_t *stokes, occa::memory &v, occa::memory &Mv)
{
  stokes->dotMultiplyKernel(stokes->Ndof,
                            v,
                            stokes->precon->invDiagA.o_v,
                            Mv);
  return;
}

static void stokesSchurComplementBlockDiagPreconditioner(stokes_t *stokes, dfloat lambda, occa::memory &v, occa::memory &Mv)
{
  occa::memory &tmp  = stokes->o_preconWorkspace[0];
  occa::memory &tmp2 = stokes->o_preconWorkspace[1];
  occa::memory &tmp3 = stokes->o_preconWorkspace[2];

  int dim = stokes->meshV->dim;

  int offset = stokes->NtotalV*sizeof(dfloat);
  
  for(int d=0;d<dim;++d){
    occa::memory  vd =  v+d*offset;
    occa::memory Mvd = Mv+d*offset;
    ellipticPreconditioner(stokes->precon->ellipticV, lambda, vd, Mvd);
  }
  
  if (stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MASSMATRIX")) {
    stokes->dotMultiplyKernel(stokes->NtotalP, stokes->precon->invMM.o_p, v+dim*offset, Mv+dim*offset);
  } else if (stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MULTIGRID")) {

    occa::memory  p =  v+dim*offset;
    occa::memory Mp = Mv+dim*offset;
    occa::memory tmpp = tmp+dim*offset;
    
    stokes->dotMultiplyKernel(stokes->NtotalP, stokes->precon->invMM.o_p, p, Mp);
    ellipticPreconditioner(stokes->precon->ellipticP, 0.0, p, tmpp);
    stokes->vecScaledAddKernel(stokes->NtotalP, lambda, tmpp, 1.0, Mp);
  }

  return;
}
