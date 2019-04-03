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
static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfoV, occa::properties& kernelInfoP);

void stokesSolveSetup(stokes_t *stokes, dfloat lambda, dfloat *eta, occa::properties &kernelInfoV, occa::properties &kernelInfoP)
{
  FILE *fp;
  char fname[BUFSIZ];
  int  Nrows, Ncols, verbose;

  verbose = stokes->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  stokes->NtotalV = stokes->meshV->Nelements*stokes->meshV->Np;
  stokes->NtotalP = stokes->meshP->Nelements*stokes->meshP->Np;
  stokes->Ndof    = stokes->meshV->dim*stokes->NtotalV + stokes->NtotalP;

  meshParallelGatherScatterSetup(stokes->meshV, stokes->NtotalV, stokes->meshV->globalIds, stokes->meshV->comm, verbose);
  meshParallelGatherScatterSetup(stokes->meshP, stokes->NtotalP, stokes->meshP->globalIds, stokes->meshP->comm, verbose);

  stokesVecAllocate(stokes, &stokes->u);
  stokesVecAllocate(stokes, &stokes->f);

  // set up work space
  stokes->NsolveWorkspace = 8; // need to check
  stokes->o_solveWorkspace = new occa::memory[stokes->NsolveWorkspace];
  for(int n=0;n<stokes->NsolveWorkspace;++n){
    stokes->o_solveWorkspace[n] = stokes->meshV->device.malloc(stokes->Ndof*sizeof(dfloat));
  }

  stokes->NpreconWorkspace = 8; // need to check
  stokes->o_preconWorkspace = new occa::memory[stokes->NpreconWorkspace];
  for(int n=0;n<stokes->NpreconWorkspace;++n){
    stokes->o_preconWorkspace[n] = stokes->meshV->device.malloc(stokes->Ndof*sizeof(dfloat));
  }
  
  stokes->eta = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  if (eta == NULL) {
    for (int i = 0; i < stokes->NtotalV; i++)
      stokes->eta[i] = 1.0;
  } else {
    for (int i = 0; i < stokes->NtotalV; i++)
      stokes->eta[i] = eta[i];
  }
  stokes->o_eta = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat), stokes->eta);

  stokes->cubEta = (dfloat*)calloc(stokes->meshV->Nelements*stokes->meshV->cubNp, sizeof(dfloat));
  for (int i = 0; i < stokes->meshV->Nelements*stokes->meshV->cubNp; i++)
    stokes->cubEta[i] = 1.0;
  stokes->o_cubEta = stokes->meshV->device.malloc(stokes->meshV->Nelements*stokes->meshV->cubNp*sizeof(dfloat), stokes->cubEta);

  if (stokes->options.compareArgs("INTEGRATION TYPE", "CUBATURE")) {
    sprintf(fname, DSTOKES "/data/stokesN%02d.dat", stokes->meshV->N);
    fp = fopen(fname, "r");
    if (!fp) {
      printf("ERROR:  Cannot open file '%s' for reading.\n", fname);
    }

    readDfloatArray(fp, "Pressure 1D cubature interpolation matrix", &stokes->cubInterpP, &Nrows, &Ncols);
    if (Nrows != stokes->meshV->cubNq) {
      printf("ERROR:  cubNq mismatch in pressure interpolation matrix (%d vs. %d).\n", Nrows, stokes->meshV->cubNq);
      exit(-1);
    } else if (Ncols != stokes->meshP->Nq) {
      printf("ERROR:  Pressure Nq mismatch in pressure interpolation matrix (%d vs. %d).\n", Ncols, stokes->meshP->Nq);
      exit(-1);
    }

    readDfloatArray(fp, "Cubature 1D differentiation matrix", &stokes->cubD, &Nrows, &Ncols);
    if (Nrows != stokes->meshV->cubNq) {
      printf("ERROR:  cubNq mismatch in differentiation matrix (%d vs. %d).\n", Nrows, stokes->meshV->cubNq);
      exit(-1);
    }

    fclose(fp);
  }

  stokes->cubInterpV = stokes->meshV->cubInterp;

  stokes->o_cubInterpV = stokes->meshV->device.malloc(stokes->meshV->Nq*stokes->meshV->cubNq*sizeof(dfloat), stokes->cubInterpV);
  stokes->o_cubInterpP = stokes->meshV->device.malloc(stokes->meshP->Nq*stokes->meshV->cubNq*sizeof(dfloat), stokes->cubInterpP);
  stokes->o_cubD = stokes->meshV->device.malloc(stokes->meshV->cubNq*stokes->meshV->cubNq*sizeof(dfloat), stokes->cubD);

  stokes->meshV->cubggeo = (dfloat*)calloc(stokes->meshV->Nelements*stokes->meshV->Nggeo*stokes->meshV->cubNp, sizeof(dfloat));
  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->cubNq; i++) {
      for (int j = 0; j < stokes->meshV->cubNq; j++) {
        dlong base;
        int n = i*stokes->meshV->cubNq + j;

        base = e*stokes->meshV->Nvgeo*stokes->meshV->cubNp + n;
        dfloat rx = stokes->meshV->cubvgeo[base + stokes->meshV->cubNp*RXID];
        dfloat ry = stokes->meshV->cubvgeo[base + stokes->meshV->cubNp*RYID];
        dfloat sx = stokes->meshV->cubvgeo[base + stokes->meshV->cubNp*SXID];
        dfloat sy = stokes->meshV->cubvgeo[base + stokes->meshV->cubNp*SYID];
        dfloat JW = stokes->meshV->cubvgeo[base + stokes->meshV->cubNp*JWID];

        base = e*stokes->meshV->Nggeo*stokes->meshV->cubNp + n;
        stokes->meshV->cubggeo[base + stokes->meshV->cubNp*G00ID] = JW*(rx*rx + ry*ry);
        stokes->meshV->cubggeo[base + stokes->meshV->cubNp*G01ID] = JW*(rx*sx + ry*sy);
        stokes->meshV->cubggeo[base + stokes->meshV->cubNp*G11ID] = JW*(sx*sx + sy*sy);
        stokes->meshV->cubggeo[base + stokes->meshV->cubNp*GWJID] = JW;
      }
    }
  }

  stokes->meshV->o_cubggeo = stokes->meshV->device.malloc(stokes->meshV->Nelements*stokes->meshV->Nggeo*stokes->meshV->cubNp*sizeof(dfloat), stokes->meshV->cubggeo);

  if (stokes->meshV->dim == 2) {
    kernelInfoV["includes"] += DSTOKES "/data/stokesBoundary2D.h";
  } else if (stokes->meshV->dim == 3) {
    kernelInfoV["includes"] += DSTOKES "/data/stokesBoundary3D.h";
  }

  stokesAllocateScratchVars(stokes);
  stokesSetupBCMask(stokes);
  stokesSetupKernels(stokes, kernelInfoV, kernelInfoP);
  stokesPreconditionerSetup(stokes, lambda, kernelInfoV, kernelInfoP);

  // set up stokes Fischer successive RHS
  int Nrhs = 0;
  stokes->options.getArgs("DIM FISCHER SUCCESSIVE RHS", Nrhs);

  stokes->fsrNrhs = 0;
  if(stokes->fsrNrhs){
    size_t NrhsEntries = (stokes->meshV->dim*stokes->NtotalV+
			  stokes->NtotalP)*stokes->fsrNrhs;

    dfloat *stokesHistoryQ = (dfloat*) calloc(NrhsEntries,sizeof(dfloat)); // can try using float here ?
    stokes->o_fsrHistoryQ = stokes->meshV->device.malloc(NrhsEntries*sizeof(dfloat),stokesHistoryQ);

    // convenient array for zeroing accumulators prior to dot products
    stokes->o_fsrZeroArray = stokes->meshV->device.malloc(stokes->fsrNrhs*sizeof(dfloat),stokesHistoryQ);
    stokes->o_fsrAlphas = stokes->meshV->device.malloc(stokes->fsrNrhs*sizeof(dfloat),stokesHistoryQ);
    
  }

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

static void stokesSetupBCMask(stokes_t *stokes)
{
  int verbose = stokes->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  // Get the BC type codes for each node, prioritizing Dirichlet over Neumann.
  stokes->mapB = (int*)calloc(stokes->NtotalV, sizeof(int));
  for (dlong e = 0; e < stokes->meshV->Nelements; e++) {
    for (int n = 0; n < stokes->meshV->Np; n++)
      stokes->mapB[e*stokes->meshV->Np + n] = INT_MAX;
    for (int f = 0; f < stokes->meshV->Nfaces; f++) {
      int bc = stokes->meshV->EToB[e*stokes->meshV->Nfaces + f];
      if (bc > 0) {
        for (int n = 0; n < stokes->meshV->Nfp; n++) {
          int fid, BCFlag;
          BCFlag = stokes->BCType[bc];
          fid = stokes->meshV->faceNodes[f*stokes->meshV->Nfp + n];
          stokes->mapB[e*stokes->meshV->Np + fid] = mymin(BCFlag, stokes->mapB[e*stokes->meshV->Np + fid]);
        }
      }
    }
  }

  ogsGatherScatter(stokes->mapB, ogsInt, ogsMin, stokes->meshV->ogs);

  // Count the number of nodes we need to mask out (Dirichlet boundary nodes).
  stokes->Nmasked = 0;
  for (dlong n = 0; n < stokes->NtotalV; n++) {
    if (stokes->mapB[n] == INT_MAX)
      stokes->mapB[n] = 0;
    else if (stokes->mapB[n] == 1)
      stokes->Nmasked++;
  }

  stokes->o_mapB = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(int), stokes->mapB);

  // Record the indices of the nodes we need to mask out.
  stokes->maskIds = (dlong*)calloc(stokes->Nmasked, sizeof(dlong));
  stokes->Nmasked = 0;
  for (dlong n = 0; n < stokes->NtotalV; n++)
    if (stokes->mapB[n] == 1)
      stokes->maskIds[stokes->Nmasked++] = n;

  if (stokes->Nmasked)
    stokes->o_maskIds = stokes->meshV->device.malloc(stokes->Nmasked*sizeof(dlong), stokes->maskIds);

  // Make the index mask.
  stokes->meshV->maskedGlobalIds = (hlong*)calloc(stokes->NtotalV, sizeof(hlong));
  memcpy(stokes->meshV->maskedGlobalIds, stokes->meshV->globalIds, stokes->NtotalV*sizeof(hlong));
  for (dlong n = 0; n < stokes->Nmasked; n++)
    stokes->meshV->maskedGlobalIds[stokes->maskIds[n]] = 0;

  // Set up the masked GS handle.
  stokes->ogs = ogsSetup(stokes->NtotalV, stokes->meshV->maskedGlobalIds, stokes->meshV->comm, verbose, stokes->meshV->device);

  return;
}

static void stokesSetupKernels(stokes_t *stokes, occa::properties &kernelInfoV, occa::properties& kernelInfoP)
{
  kernelInfoV["defines/p_blockSize"] = STOKES_REDUCTION_BLOCK_SIZE;
  kernelInfoV["defines/p_DIM"]   = stokes->meshV->dim;
  kernelInfoV["defines/p_NpV"]   = stokes->meshV->Np;
  kernelInfoV["defines/p_NqV"]   = stokes->meshV->Nq;
  kernelInfoV["defines/p_NpP"]   = stokes->meshP->Np;
  kernelInfoV["defines/p_NqP"]   = stokes->meshP->Nq;

  int fsrNrhs = 0;
  stokes->options.getArgs("DIM FISCHER SUCCESSIVE RHS", fsrNrhs);
  kernelInfoV["defines/p_fsrNrhs"] = fsrNrhs;
  
  stokes->meshV->maskKernel          = stokes->meshV->device.buildKernel(DHOLMES "/okl/mask.okl", "mask", kernelInfoV);

  stokes->dotMultiplyKernel          = stokes->meshV->device.buildKernel(DHOLMES "/okl/dotMultiply.okl", "dotMultiply", kernelInfoV);
  stokes->vecScaleKernel             = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecScale.okl", "stokesVecScale", kernelInfoV);
  stokes->vecScaledAddKernel         = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecScaledAdd.okl", "stokesVecScaledAdd", kernelInfoV);
  stokes->vecZeroKernel              = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesVecZero.okl", "stokesVecZero", kernelInfoV);
  stokes->weightedInnerProductKernel = stokes->meshV->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl", "weightedInnerProduct2", kernelInfoV);

  stokes->globalWeightedInnerProductKernel
    = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesGlobalWeightedInnerProduct.okl", "stokesGlobalWeightedInnerProduct", kernelInfoV);


  stokes->multipleGlobalWeightedInnerProductsKernel
    = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesMultipleGlobalWeightedInnerProducts.okl", "stokesMultipleGlobalWeightedInnerProducts", kernelInfoV);


  stokes->fsrStartKernel
    = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesFSRStart.okl", "stokesFSRStart", kernelInfoV);


  stokes->fsrUpdateKernel
    = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesFSRUpdate.okl", "stokesFSRUpdate", kernelInfoV);
  
  /* TODO:  Replace this with parametrized filenames. */
  if ((stokes->meshV->dim == 2) && (stokes->elementType == QUADRILATERALS)) {
    stokes->lowerPressureKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesLowerPressureQuad2D.okl", "stokesLowerPressureQuad2D", kernelInfoV);
    stokes->raisePressureKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesRaisePressureQuad2D.okl", "stokesRaisePressureQuad2D", kernelInfoV);

    if (stokes->options.compareArgs("INTEGRATION TYPE", "GLL")) {
      stokes->divergenceKernel     = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesDivergenceQuad2D.okl", "stokesDivergenceQuad2D", kernelInfoV);
      stokes->gradientKernel       = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesGradientQuad2D.okl", "stokesGradientQuad2D", kernelInfoV);
      stokes->stiffnessKernel      = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesStiffnessQuad2D.okl", "stokesStiffnessQuad2D", kernelInfoV);
      stokes->stokesOperatorKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesOperatorQuad2D.okl", "stokesOperatorQuad2D", kernelInfoV);
    } else if (stokes->options.compareArgs("INTEGRATION TYPE", "CUBATURE")) {
      stokes->divergenceKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesDivergenceQuad2D.okl", "stokesDivergenceCubatureQuad2D", kernelInfoV);
      stokes->gradientKernel   = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesGradientQuad2D.okl", "stokesGradientCubatureQuad2D", kernelInfoV);
      stokes->stiffnessKernel  = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesStiffnessQuad2D.okl", "stokesStiffnessCubatureQuad2D", kernelInfoV);
    } else {
      printf("ERROR:  Invalid value for option [INTEGRATION TYPE].\n");
      exit(-1);
    }

  } else if ((stokes->meshV->dim == 3) && (stokes->elementType == HEXAHEDRA)) {
    stokes->lowerPressureKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesLowerPressureHex3D.okl", "stokesLowerPressureHex3D", kernelInfoV);
    stokes->raisePressureKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesRaisePressureHex3D.okl", "stokesRaisePressureHex3D", kernelInfoV);

    if (stokes->options.compareArgs("INTEGRATION TYPE", "GLL")) {
      stokes->divergenceKernel     = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesDivergenceHex3D.okl", "stokesDivergenceHex3D", kernelInfoV);
      stokes->gradientKernel       = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesGradientHex3D.okl", "stokesGradientHex3D", kernelInfoV);
      stokes->stiffnessKernel      = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesStiffnessHex3D.okl", "stokesStiffnessHex3D", kernelInfoV);
      stokes->stokesOperatorKernel = stokes->meshV->device.buildKernel(DSTOKES "/okl/stokesOperatorHex3D.okl", "stokesOperatorHex3D", kernelInfoV);
    } else if (stokes->options.compareArgs("INTEGRATION TYPE", "CUBATURE")) {
      printf("ERROR:  Not implemented.\n");
      MPI_Finalize();
      exit(-1);
    }
  }
  
  return;
}
