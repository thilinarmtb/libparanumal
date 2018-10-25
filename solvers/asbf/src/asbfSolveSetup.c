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

#include "asbf.h"

void asbfSolveSetup(asbf_t *asbf, dfloat lambda, occa::properties &kernelInfo)
{
  mesh_t *mesh = asbf->mesh;
  setupAide options = asbf->options;

  // TODO:  Compute these matrices on-the-fly.
  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/solvers/asbf/data/asbfN%02d.dat", asbf->Nmodes);

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  int Nrows, Ncols;

  readDfloatArray(fp, "ASBF EIGENVALUES",
      &(asbf->eigenvalues),&(asbf->Nmodes), &(Ncols));
  readDfloatArray(fp, "ASBF QUADRATURE VANDERMONDE",
      &(asbf->Bquad),&(asbf->Nquad), &(asbf->Nmodes));
  readDfloatArray(fp, "ASBF QUADRATURE NODES",
      &(asbf->Rquad),&(asbf->Nquad), &Ncols);
  readDfloatArray(fp, "ASBF QUADRATURE WEIGHTS",
      &(asbf->Wquad),&(asbf->Nquad), &Ncols);
  readDfloatArray(fp, "ASBF GLL VANDERMONDE",
      &(asbf->Bgll),&(asbf->Ngll), &(asbf->Nmodes));
  readDfloatArray(fp, "ASBF GLL NODES",
      &(asbf->Rgll),&(asbf->Ngll), &Ncols);
  readDfloatArray(fp, "ASBF PLOT VANDERMONDE",
      &(asbf->Bplot),&(asbf->Nplot), &(asbf->Nmodes));
  readDfloatArray(fp, "ASBF PLOT NODES",
      &(asbf->Rplot),&(asbf->Nplot), &Ncols);
  readDfloatArray(fp, "ASBF QUADRATURE DERIVATIVE VANDERMONDE",
      &(asbf->DBquad),&(asbf->Nquad), &(asbf->Nmodes));
  readDfloatArray(fp, "ASBF GLL DERIVATIVE VANDERMONDE",
      &(asbf->DBgll),&(asbf->Ngll), &(asbf->Nmodes));
  readDfloatArray(fp, "ASBF PLOT DERIVATIVE VANDERMONDE",
      &(asbf->DBplot),&(asbf->Nplot), &(asbf->Nmodes));

  fclose(fp);

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
  asbf->Ntotal = Nlocal + Nhalo;

  asbf->r3D = (dfloat*) calloc(asbf->Nmodes*asbf->Ntotal, sizeof(dfloat));
  asbf->q3D = (dfloat*) calloc(asbf->Nmodes*asbf->Ntotal, sizeof(dfloat));
  asbf->f   = (dfloat*) calloc(asbf->Nquad, sizeof(dfloat));
  asbf->r   = (dfloat*) calloc(asbf->Ntotal, sizeof(dfloat));
  asbf->x   = (dfloat*) calloc(asbf->Ntotal, sizeof(dfloat));

  kernelInfo["defines/p_asbfNmodes"] = asbf->Nmodes;
  kernelInfo["defines/p_asbfNquad"] = asbf->Nquad;
  kernelInfo["defines/p_asbfNgll"] = asbf->Ngll;

  occa::properties kernelInfoP  = kernelInfo;

  asbf->o_r   = mesh->device.malloc(asbf->Ntotal*sizeof(dfloat), asbf->r);
  asbf->o_x   = mesh->device.malloc(asbf->Ntotal*sizeof(dfloat), asbf->x);

  asbf->pOptions = options;

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int pBCType[7] = {0,1,1,2,1,1,1}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances
  asbf->pTOL = 1E-8;

  asbf->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  asbf->pSolver->mesh = mesh;
  asbf->pSolver->options = asbf->pOptions;
  asbf->pSolver->dim = asbf->dim;
  asbf->pSolver->elementType = asbf->elementType;
  asbf->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(asbf->pSolver->BCType,pBCType,7*sizeof(int));

  ellipticSolveSetup(asbf->pSolver, asbf->lambda, kernelInfoP);

  asbf->meshSEM = NULL;

  if(asbf->elementType==QUADRILATERALS){
    asbfExtrudeSphere(asbf);
  }


  // OKL kernels specific to asbf
  asbf->asbfReconstructKernel =
    mesh->device.buildKernel(DASBF "/okl/asbfReconstructHex3D.okl",
        "asbfReconstructHex3D",
        kernelInfo);
}
