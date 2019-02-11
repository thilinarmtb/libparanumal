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

static void stokesAllocateScratchVars(stokes_t *stokes);
static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfoV, occa::properties& kernelInfoP);

void stokesSolveSetup(stokes_t *stokes, dfloat *eta, occa::properties &kernelInfoV, occa::properties &kernelInfoP)
{
  int verbose = stokes->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  stokes->NtotalV = stokes->meshV->Nelements*stokes->meshV->Np;
  stokes->NtotalP = stokes->meshP->Nelements*stokes->meshP->Np;
  stokes->Ndof    = stokes->meshV->dim*stokes->NtotalV + stokes->NtotalP;

  meshParallelGatherScatterSetup(stokes->meshV, stokes->NtotalV, stokes->meshV->globalIds, stokes->meshV->comm, verbose);
  meshParallelGatherScatterSetup(stokes->meshP, stokes->NtotalP, stokes->meshP->globalIds, stokes->meshP->comm, verbose);

  stokesVecAllocate(stokes, &stokes->u);
  stokesVecAllocate(stokes, &stokes->f);

  stokes->eta = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  if (eta == NULL) {
    for (int i = 0; i < stokes->NtotalV; i++)
      stokes->eta[i] = 1.0;
  } else {
    for (int i = 0; i < stokes->NtotalV; i++)
      stokes->eta[i] = eta[i];
  }

  stokes->o_eta = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat), stokes->eta);

  stokesAllocateScratchVars(stokes);
  stokesSetupKernels(stokes, kernelInfoV, kernelInfoP);
  stokesPreconditionerSetup(stokes);

  return;
}

static void stokesAllocateScratchVars(stokes_t *stokes)
{
  stokes->NblockV = mymax(1, (stokes->NtotalV + STOKES_REDUCTION_BLOCK_SIZE - 1)/STOKES_REDUCTION_BLOCK_SIZE);
  stokes->workV = (dfloat*)calloc(stokes->NblockV, sizeof(dfloat));
  stokes->o_workV = stokes->meshV->device.malloc(stokes->NblockV*sizeof(dfloat), stokes->workV);

  stokes->NblockP = mymax(1, (stokes->NtotalP + STOKES_REDUCTION_BLOCK_SIZE - 1)/STOKES_REDUCTION_BLOCK_SIZE);
  stokes->workP = (dfloat*)calloc(stokes->NblockP, sizeof(dfloat));
  stokes->o_workP = stokes->meshV->device.malloc(stokes->NblockP*sizeof(dfloat), stokes->workP);

  return;
}

static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfoV, occa::properties& kernelInfoP)
{
  kernelInfoV["defines/p_blockSize"] = STOKES_REDUCTION_BLOCK_SIZE;
  kernelInfoV["defines/p_NpV"] = stokes->meshV->Np;
  kernelInfoV["defines/p_NqV"] = stokes->meshV->Nq;
  kernelInfoV["defines/p_NpP"] = stokes->meshP->Np;
  kernelInfoV["defines/p_NqP"] = stokes->meshP->Nq;

  stokes->divergenceKernel           = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesDivergenceQuad2D.okl", "stokesDivergenceQuad2D", kernelInfoV);
  stokes->dotMultiplyKernel          = stokes->meshV->device.buildKernel(DHOLMES "/okl/dotMultiply.okl", "dotMultiply", kernelInfoV);
  stokes->gradientKernel             = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesGradientQuad2D.okl", "stokesGradientQuad2D", kernelInfoV);
  stokes->lowerPressureKernel        = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesLowerPressureQuad2D.okl", "stokesLowerPressureQuad2D", kernelInfoV);
  stokes->raisePressureKernel        = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesRaisePressureQuad2D.okl", "stokesRaisePressureQuad2D", kernelInfoV);
  stokes->stiffnessKernel            = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesStiffnessQuad2D.okl", "stokesStiffnessQuad2D", kernelInfoV);
  stokes->vecScaleKernel             = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecScale.okl", "stokesVecScale", kernelInfoV);
  stokes->vecScaledAddKernel         = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecScaledAdd.okl", "stokesVecScaledAdd", kernelInfoV);
  stokes->vecZeroKernel              = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecZero.okl", "stokesVecZero", kernelInfoV);
  stokes->weightedInnerProductKernel = stokes->meshV->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl", "weightedInnerProduct2", kernelInfoV);

  return;
}
