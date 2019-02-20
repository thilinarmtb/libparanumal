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

#include <limits.h>
#include "stokes.h"

static void stokesAllocateScratchVars(stokes_t *stokes);
static void stokesSetupBCMask(stokes_t *stokes);
static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfo);

void stokesSolveSetup(stokes_t *stokes, dfloat *eta, occa::properties &kernelInfo)
{
  FILE *fp;
  char fname[BUFSIZ];
  int  Nrows, Ncols, verbose;

  verbose = stokes->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  stokes->Ntotal = stokes->mesh->Nelements*stokes->mesh->Np;
  stokes->Ndof   = (stokes->mesh->dim + 1)*stokes->Ntotal;

  meshParallelGatherScatterSetup(stokes->mesh, stokes->Ntotal, stokes->mesh->globalIds, stokes->mesh->comm, verbose);

  stokesVecAllocate(stokes, &stokes->u);
  stokesVecAllocate(stokes, &stokes->f);

  stokes->eta = (dfloat*)calloc(stokes->Ntotal, sizeof(dfloat));
  if (eta == NULL) {
    for (int i = 0; i < stokes->Ntotal; i++)
      stokes->eta[i] = 1.0;
  } else {
    for (int i = 0; i < stokes->Ntotal; i++)
      stokes->eta[i] = eta[i];
  }

  stokes->o_eta = stokes->mesh->device.malloc(stokes->Ntotal*sizeof(dfloat), stokes->eta);

  if (stokes->mesh->dim == 2) {
    kernelInfo["includes"] += DSTOKES "/data/stokesBoundary2D.h";
  } else if (stokes->mesh->dim == 3) {
    printf("ERROR:  Not implemented.\n");
    exit(-1);
  }

  sprintf(fname, DSTOKES "/data/stokesN%02d.dat", stokes->mesh->N);
  fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR:  Cannot open file '%s' for reading.\n", fname);
  }

  readDfloatArray(fp, "Pressure projection rank-1 update - Interpolatory (u)", &stokes->uP, &Nrows, &Ncols);
  readDfloatArray(fp, "Pressure projection rank-1 update - Interpolatory (v)", &stokes->vP, &Nrows, &Ncols);

  fclose(fp);

  stokes->o_uP = stokes->mesh->device.malloc(stokes->mesh->Nq*sizeof(dfloat), stokes->uP);
  stokes->o_vP = stokes->mesh->device.malloc(stokes->mesh->Nq*sizeof(dfloat), stokes->vP);

  stokesAllocateScratchVars(stokes);
  stokesSetupBCMask(stokes);
  stokesSetupKernels(stokes, kernelInfo);
  stokesPreconditionerSetup(stokes);

  return;
}

static void stokesAllocateScratchVars(stokes_t *stokes)
{
  stokes->Nblock = mymax(1, (stokes->Ntotal + STOKES_REDUCTION_BLOCK_SIZE - 1)/STOKES_REDUCTION_BLOCK_SIZE);
  stokes->block = (dfloat*)calloc(stokes->Nblock, sizeof(dfloat));
  stokes->o_block = stokes->mesh->device.malloc(stokes->Nblock*sizeof(dfloat), stokes->block);

  return;
}

static void stokesSetupBCMask(stokes_t *stokes)
{
  int verbose = stokes->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  // Get the BC type codes for each node, prioritizing Dirichlet over Neumann.
  stokes->mapB = (int*)calloc(stokes->Ntotal, sizeof(int));
  for (dlong e = 0; e < stokes->mesh->Nelements; e++) {
    for (int n = 0; n < stokes->mesh->Np; n++)
      stokes->mapB[e*stokes->mesh->Np + n] = INT_MAX;
    for (int f = 0; f < stokes->mesh->Nfaces; f++) {
      int bc = stokes->mesh->EToB[e*stokes->mesh->Nfaces + f];
      if (bc > 0) {
        for (int n = 0; n < stokes->mesh->Nfp; n++) {
          int fid, BCFlag;
          BCFlag = stokes->BCType[bc];
          fid = stokes->mesh->faceNodes[f*stokes->mesh->Nfp + n];
          stokes->mapB[e*stokes->mesh->Np + fid] = mymin(BCFlag, stokes->mapB[e*stokes->mesh->Np + fid]);
        }
      }
    }
  }

  ogsGatherScatter(stokes->mapB, ogsInt, ogsMin, stokes->mesh->ogs);

  // Count the number of nodes we need to mask out (Dirichlet boundary nodes).
  stokes->Nmasked = 0;
  for (dlong n = 0; n < stokes->Ntotal; n++) {
    if (stokes->mapB[n] == INT_MAX)
      stokes->mapB[n] = 0;
    else if (stokes->mapB[n] == 1)
      stokes->Nmasked++;
  }

  stokes->o_mapB = stokes->mesh->device.malloc(stokes->Ntotal*sizeof(int), stokes->mapB);

  // Record the indices of the nodes we need to mask out.
  stokes->maskIds = (dlong*)calloc(stokes->Nmasked, sizeof(dlong));
  stokes->Nmasked = 0;
  for (dlong n = 0; n < stokes->Ntotal; n++)
    if (stokes->mapB[n] == 1)
      stokes->maskIds[stokes->Nmasked++] = n;

  if (stokes->Nmasked)
    stokes->o_maskIds = stokes->mesh->device.malloc(stokes->Nmasked*sizeof(dlong), stokes->maskIds);

  // Make the index mask.
  stokes->mesh->maskedGlobalIds = (hlong*)calloc(stokes->Ntotal, sizeof(hlong));
  memcpy(stokes->mesh->maskedGlobalIds, stokes->mesh->globalIds, stokes->Ntotal*sizeof(hlong));
  for (dlong n = 0; n < stokes->Nmasked; n++)
    stokes->mesh->maskedGlobalIds[stokes->maskIds[n]] = 0;

  // Set up the masked GS handle.
  stokes->ogs = ogsSetup(stokes->Ntotal, stokes->mesh->maskedGlobalIds, stokes->mesh->comm, verbose, stokes->mesh->device);

  return;
}

static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfo)
{
  kernelInfo["defines/p_blockSize"] = STOKES_REDUCTION_BLOCK_SIZE;

  stokes->mesh->maskKernel          = stokes->mesh->device.buildKernel(DHOLMES "/okl/mask.okl", "mask", kernelInfo);

  stokes->dotMultiplyKernel          = stokes->mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl", "dotMultiply", kernelInfo);
  stokes->vecScaleKernel             = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesVecScale.okl", "stokesVecScale", kernelInfo);
  stokes->vecScaledAddKernel         = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesVecScaledAdd.okl", "stokesVecScaledAdd", kernelInfo);
  stokes->vecZeroKernel              = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesVecZero.okl", "stokesVecZero", kernelInfo);
  stokes->weightedInnerProductKernel = stokes->mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl", "weightedInnerProduct2", kernelInfo);
  stokes->updateMINRESKernel         = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesUpdateMINRES.okl", "stokesUpdateMINRES", kernelInfo);
  
  /* TODO:  Replace this with parametrized filenames. */
  if ((stokes->mesh->dim == 2) && (stokes->elementType == QUADRILATERALS)) {
    stokes->divergenceKernel        = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesDivergenceQuad2D.okl", "stokesDivergenceQuad2D", kernelInfo);
    stokes->gradientKernel          = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesGradientQuad2D.okl", "stokesGradientQuad2D", kernelInfo);
    stokes->rankOneProjectionKernel = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesRankOneProjectionQuad2D.okl", "stokesRankOneProjectionQuad2D", kernelInfo);
    stokes->stiffnessKernel         = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesStiffnessQuad2D.okl", "stokesStiffnessQuad2D", kernelInfo);
  } else if ((stokes->mesh->dim == 3) && (stokes->elementType == HEXAHEDRA)) {
    stokes->divergenceKernel        = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesDivergenceHex3D.okl", "stokesDivergenceHex3D", kernelInfo);
    stokes->gradientKernel          = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesGradientHex3D.okl", "stokesGradientHex3D", kernelInfo);
    stokes->rankOneProjectionKernel = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesRankOneProjectionHex3D.okl", "stokesRankOneProjectionHex3D", kernelInfo);
    stokes->stiffnessKernel         = stokes->mesh->device.buildKernel(DSTOKES "/okl/stokesStiffnessHex3D.okl", "stokesStiffnessHex3D", kernelInfo);
  }

  return;
}
