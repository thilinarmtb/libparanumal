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

extern "C"
{
  void dsygv_(int    *ITYPE,
              char   *JOBZ,
              char   *UPLO,
              int    *N,
              double *A,
              int    *LDA,
              double *B,
              int    *LDB,
              double *W,
              double *WORK,
              int    *LWORK,
              int    *INFO);
}

static void asbfScaleNodesAndWeights(asbf_t *asbf);
static void asbfAssembleRadialVandermondeMatrices(asbf_t *asbf, dfloat lambda, FILE *fp);

void asbfSolveSetup(asbf_t *asbf, dfloat lambda, occa::properties &kernelInfo)
{
  mesh_t *mesh = asbf->mesh;
  setupAide options = asbf->options;
  asbf->lambda = lambda;

  options.getArgs("OUTER RADIUS", asbf->R);
  if (asbf->R <= 1) {
    printf("ERROR:  OUTER RADIUS must be greater than 1.\n");
    exit(-1);
  }

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/solvers/asbf/data/asbfN%02d.dat", asbf->Nmodes);

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  int Nrows, Ncols;

  readDfloatArray(fp, "ASBF QUADRATURE NODES",
      &(asbf->Rquad),&(asbf->Nquad), &Ncols);
  readDfloatArray(fp, "ASBF QUADRATURE WEIGHTS",
      &(asbf->Wquad),&(asbf->Nquad), &Ncols);
  readDfloatArray(fp, "ASBF GLL NODES",
      &(asbf->Rgll),&(asbf->Ngll), &Ncols);
  readDfloatArray(fp, "ASBF PLOT NODES",
      &(asbf->Rplot),&(asbf->Nplot), &Ncols);

  asbfScaleNodesAndWeights(asbf);
  asbfAssembleRadialVandermondeMatrices(asbf, lambda, fp);

  // Put the r^2 factor into the quadrature weights.
  //
  // TODO:  Fix other codes so this can be removed.
  for (int i = 0; i < asbf->Nquad; i++)
    asbf->Wquad[i] *= pow(asbf->Rquad[i], 2);

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

  asbf->elliptic = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  asbf->elliptic->mesh = mesh;
  asbf->elliptic->options = asbf->options;
  asbf->elliptic->dim = asbf->dim;
  asbf->elliptic->elementType = asbf->elementType;
  asbf->elliptic->BCType = (int*) calloc(7,sizeof(int));
  memcpy(asbf->elliptic->BCType,pBCType,7*sizeof(int));

  ellipticSolveSetup(asbf->elliptic, asbf->lambda, kernelInfoP);

  asbf->meshSEM = NULL;

  // TODO:  asbfExtrudeSphere needs to be rewritten to remove the assumption
  // that asbf->Ngll == mesh->Nq.  We don't use asbf->meshSEM for anything
  // anyway, so disable this for now.
  //
  // if(asbf->elementType==QUADRILATERALS){
  //   asbfExtrudeSphere(asbf);
  // }

  // OKL kernels specific to asbf
  asbf->asbfReconstructKernel =
    mesh->device.buildKernel(DASBF "/okl/asbfReconstructHex3D.okl",
        "asbfReconstructHex3D",
        kernelInfo);
}


/*****************************************************************************/

// Scale node sets and quadrature weights from [-1, 1] to [1, R].
static void asbfScaleNodesAndWeights(asbf_t *asbf)
{
  for (int i = 0; i < asbf->Nquad; i++) {
    asbf->Rquad[i] = (1 + asbf->Rquad[i])*((asbf->R - 1)/2) + 1;
    asbf->Wquad[i] = asbf->Wquad[i]*(asbf->R - 1)/2;
  }

  for (int i = 0; i < asbf->Ngll; i++)
    asbf->Rgll[i] = (1 + asbf->Rgll[i])*((asbf->R - 1)/2) + 1;
  for (int i = 0; i < asbf->Nplot; i++)
    asbf->Rplot[i] = (1 + asbf->Rplot[i])*((asbf->R - 1)/2) + 1;
}

// Compute the Vandermonde matrices for the radial basis functions and their
// derivatives over the various node sets.  This must be called AFTER the node
// sets have been read from the setup file and scaled into [1, R].
//
// TODO:  Neumann boundary conditions.
static void asbfAssembleRadialVandermondeMatrices(asbf_t *asbf, dfloat lambda, FILE *fp)
{
  dfloat *Tquad, *Tgll, *Tplot;    // Initial basis Vandermonde matrices.
  dfloat *DTquad, *DTgll, *DTplot; // Initial basis derivative Vandermondes.
  dfloat *A, *B;                   // Matrices for GEVP.
  dfloat C1, C2;                   // Scaling constants (e.g., Jacobian).
  int Nrows, Ncols;                // Needed by readDfloatArray().
  int Nmodes, Nquad, Ngll, Nplot;  // Local copies of variables from asbf.

  // Variables needed for dsygv_()---see LAPACK documentation.
  int ITYPE, N, LDA, LDB, LWORK, INFO;
  char JOBZ, UPLO;
  dfloat *W, *WORK;

  Nmodes = asbf->Nmodes;
  Nquad  = asbf->Nquad;
  Ngll   = asbf->Ngll;
  Nplot  = asbf->Nplot;

  readDfloatArray(fp, "ASBF INITIAL BASIS QUADRATURE VANDERMONDE",
    &Tquad, &Nrows, &Ncols);
  readDfloatArray(fp, "ASBF INITIAL BASIS QUADRATURE DERIVATIVE VANDERMONDE",
    &DTquad, &Nrows, &Ncols);
  readDfloatArray(fp, "ASBF INITIAL BASIS GLL VANDERMONDE",
    &Tgll, &Nrows, &Ncols);
  readDfloatArray(fp, "ASBF INITIAL BASIS GLL DERIVATIVE VANDERMONDE",
    &DTgll, &Nrows, &Ncols);
  readDfloatArray(fp, "ASBF INITIAL BASIS PLOT VANDERMONDE",
    &Tplot, &Nrows, &Ncols);
  readDfloatArray(fp, "ASBF INITIAL BASIS PLOT DERIVATIVE VANDERMONDE",
    &DTplot, &Nrows, &Ncols);

  // Scale initial basis matrices from [-1, 1] to [1, R].
  //
  // TODO:  We need to scale DTxxx, but do we need to scale Txxx too?
  // TODO:  The C1 scaling only applies for Dirichlet conditions.
  C1 = pow((asbf->R - 1)/2, 2);
  C2 = (asbf->R - 1)/2;

  for (int i = 0; i < Nquad*Nmodes; i++) {
      Tquad[i] *= C1;
      DTquad[i] *= C2;
  }

  for (int i = 0; i < Ngll*Nmodes; i++) {
      Tgll[i]  *= C1;
      DTgll[i]  *= C2;
  }

  for (int i = 0; i < Nplot*Nmodes; i++) {
      Tplot[i] *= C1;
      DTplot[i] *= C2;
  }

  // Set up and solve the eigenvalue problem.
  A = (dfloat*) calloc(Nmodes*Nmodes, sizeof(dfloat)); 
  B = (dfloat*) calloc(Nmodes*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Nmodes; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nquad; k++) {
        dfloat a1, a2, b;
        a1 = asbf->Wquad[k]*pow(asbf->Rquad[k], 2)*DTquad[i + k*Nmodes]*DTquad[j + k*Nmodes];
        a2 = asbf->lambda*asbf->Wquad[k]*pow(asbf->Rquad[k], 2)*Tquad[i + k*Nmodes]*Tquad[j + k*Nmodes];
        b = asbf->Wquad[k]*Tquad[i + k*Nmodes]*Tquad[j + k*Nmodes];
        A[i*Nmodes + j] += a1 + a2;
        B[i*Nmodes + j] += b;
      }
    }
  }

  ITYPE = 1;        // Solve Ax = cBx
  JOBZ = 'V';       // Want both eigenvalues and eigenvectors.
  UPLO = 'U';       // Assume upper triangular storage for symmetric matrices.
  N = Nmodes;       // Matrix dimension.
  LDA = Nmodes;     // Leading dimension of A.
  LDB = Nmodes;     // Leading dimension of B.
  LWORK = 4*Nmodes; // Workspace size.  (TODO:  ilaenv() for optimal value?)

  W = (dfloat*) calloc(N, sizeof(dfloat));        // Eigenvalues.
  WORK = (dfloat*) calloc(LWORK, sizeof(dfloat)); // LAPACK workspace.

  dsygv_(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);

  if (INFO != 0)
    printf("Error in DSYGV (INFO = %d).\n", INFO);

  // Assign eigenvalues.  (TODO:  Sort them and the modes first?)
  asbf->eigenvalues = (dfloat*) calloc(Nmodes, sizeof(dfloat));
  memcpy(asbf->eigenvalues, W, Nmodes*sizeof(dfloat));

  // Assemble Vandermonde matrices.
  //
  // NB:  A contains the eigenvectors, stored in column-major order.
  asbf->Bquad = (dfloat*) calloc(Nquad*Nmodes, sizeof(dfloat));
  asbf->DBquad = (dfloat*) calloc(Nquad*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Nquad; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        asbf->Bquad[i*Nmodes + j]  += Tquad[i*Nmodes + k]*A[j*Nmodes + k];
        asbf->DBquad[i*Nmodes + j] += DTquad[i*Nmodes + k]*A[j*Nmodes + k];
      }
    }
  }

  asbf->Bgll = (dfloat*) calloc(Ngll*Nmodes, sizeof(dfloat));
  asbf->DBgll = (dfloat*) calloc(Ngll*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Ngll; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        asbf->Bgll[i*Nmodes + j]  += Tgll[i*Nmodes + k]*A[j*Nmodes + k];
        asbf->DBgll[i*Nmodes + j] += DTgll[i*Nmodes + k]*A[j*Nmodes + k];
      }
    }
  }

  asbf->Bplot = (dfloat*) calloc(Nplot*Nmodes, sizeof(dfloat));
  asbf->DBplot = (dfloat*) calloc(Nplot*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Nplot; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        asbf->Bplot[i*Nmodes + j]   += Tplot[i*Nmodes + k]*A[j*Nmodes + k];
        asbf->DBplot[i*Nmodes + j]  += DTplot[i*Nmodes + k]*A[j*Nmodes + k];
      }
    }
  }

  // We're done.  Clean up.
  free(A);
  free(B);
  free(W);
  free(WORK);
  free(Tquad);
  free(Tgll);
  free(Tplot);
  free(DTquad);
  free(DTgll);
  free(DTplot);
}
