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

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/solvers/asbf/data/asbfN%02d.dat", asbf->Nmodes);

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  int Nrows, Ncols;

  // TODO:  Read outer radius from setup file.
  asbf->R = 1.5;

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

  // TODO:  Fix this off-by-one error.
  asbf->Nmodes += 1;

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

  if(asbf->elementType==QUADRILATERALS){
    asbfExtrudeSphere(asbf);
  }

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
  int Ndeg;                        // Expansion degree (no. of modes + 1).
  int Nrows, Ncols;                // Needed by readDfloatArray().

  // Variables needed for dsygv_()---see LAPACK documentation.
  int ITYPE, N, LDA, LDB, LWORK, INFO;
  char JOBZ, UPLO;
  dfloat *W, *WORK;

  // TODO:  Fix the off-by-one error so that this parameter is not necessary.
  Ndeg = asbf->Nmodes + 1;

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

  for (int i = 0; i < asbf->Nquad*Ndeg; i++) {
      Tquad[i] *= C1;
      DTquad[i] *= C2;
  }

  for (int i = 0; i < asbf->Ngll*Ndeg; i++) {
      Tgll[i]  *= C1;
      DTgll[i]  *= C2;
  }

  for (int i = 0; i < asbf->Nplot*Ndeg; i++) {
      Tplot[i] *= C1;
      DTplot[i] *= C2;
  }

  // Set up and solve the eigenvalue problem.
  A = (dfloat*) calloc(Ndeg*Ndeg, sizeof(dfloat)); 
  B = (dfloat*) calloc(Ndeg*Ndeg, sizeof(dfloat));
  for (int i = 0; i < Ndeg; i++) {
    for (int j = 0; j < Ndeg; j++) {
      for (int k = 0; k < asbf->Nquad; k++) {
        dfloat a1, a2, b;
        a1 = asbf->Wquad[k]*pow(asbf->Rquad[k], 2)*DTquad[i + k*Ndeg]*DTquad[j + k*Ndeg];
        a2 = asbf->lambda*asbf->Wquad[k]*pow(asbf->Rquad[k], 2)*Tquad[i + k*Ndeg]*Tquad[j + k*Ndeg];
        b = asbf->Wquad[k]*Tquad[i + k*Ndeg]*Tquad[j + k*Ndeg];
        A[i*Ndeg + j] += a1 + a2;
        B[i*Ndeg + j] += b;
      }
    }
  }

  ITYPE = 1;      // Solve Ax = cBx
  JOBZ = 'V';     // Want both eigenvalues and eigenvectors.
  UPLO = 'U';     // Assume upper triangular storage for symmetric matrices.
  N = Ndeg;       // Matrix dimension.
  LDA = Ndeg;     // Leading dimension of A.
  LDB = Ndeg;     // Leading dimension of B.
  LWORK = 4*Ndeg; // Workspace size.  (TODO:  ilaenv() for optimal value?)

  W = (dfloat*) calloc(N, sizeof(dfloat));        // Eigenvalues.
  WORK = (dfloat*) calloc(LWORK, sizeof(dfloat)); // LAPACK workspace.

  dsygv_(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);

  if (INFO != 0)
    printf("Error in DSYGV (INFO = %d).\n", INFO);

  // Assign eigenvalues.  (TODO:  Sort them and the modes first?)
  asbf->eigenvalues = (dfloat*) calloc(Ndeg, sizeof(dfloat));
  memcpy(asbf->eigenvalues, W, Ndeg*sizeof(dfloat));

  // Assemble Vandermonde matrices.
  //
  // NB:  A contains the eigenvectors, stored in column-major order.
  asbf->Bquad = (dfloat*) calloc(asbf->Nquad*Ndeg, sizeof(dfloat));
  asbf->DBquad = (dfloat*) calloc(asbf->Nquad*Ndeg, sizeof(dfloat));
  for (int i = 0; i < asbf->Nquad; i++) {
    for (int j = 0; j < Ndeg; j++) {
      for (int k = 0; k < Ndeg; k++) {
        asbf->Bquad[i*Ndeg + j]  += Tquad[i*Ndeg + k]*A[j*Ndeg + k];
        asbf->DBquad[i*Ndeg + j] += DTquad[i*Ndeg + k]*A[j*Ndeg + k];
      }
    }
  }

  asbf->Bgll = (dfloat*) calloc(asbf->Ngll*Ndeg, sizeof(dfloat));
  asbf->DBgll = (dfloat*) calloc(asbf->Ngll*Ndeg, sizeof(dfloat));
  for (int i = 0; i < asbf->Ngll; i++) {
    for (int j = 0; j < Ndeg; j++) {
      for (int k = 0; k < Ndeg; k++) {
        asbf->Bgll[i*Ndeg + j]  += Tgll[i*Ndeg + k]*A[j*Ndeg + k];
        asbf->DBgll[i*Ndeg + j] += DTgll[i*Ndeg + k]*A[j*Ndeg + k];
      }
    }
  }

  asbf->Bplot = (dfloat*) calloc(asbf->Nplot*Ndeg, sizeof(dfloat));
  asbf->DBplot = (dfloat*) calloc(asbf->Nplot*Ndeg, sizeof(dfloat));
  for (int i = 0; i < asbf->Nplot; i++) {
    for (int j = 0; j < Ndeg; j++) {
      for (int k = 0; k < Ndeg; k++) {
        asbf->Bplot[i*Ndeg + j]   += Tplot[i*Ndeg + k]*A[j*Ndeg + k];
        asbf->DBplot[i*Ndeg + j]  += DTplot[i*Ndeg + k]*A[j*Ndeg + k];
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
