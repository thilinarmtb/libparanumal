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

#include "shell.h"

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

static void shellLoadRadialBasisTrue(shell_t *shell);
static void shellLoadRadialBasisGlobalDiscrete(shell_t *shell);
static void shellLoadRadialBasisPiecewiseDiscrete(shell_t *shell);

void shellSolveSetup(shell_t *shell, dfloat lambda, occa::properties &kernelInfo)
{
  mesh_t *mesh = shell->mesh;
  setupAide options = shell->options;
  shell->lambda = lambda;

  options.getArgs("OUTER RADIUS", shell->R);
  if (shell->R <= 1) {
    printf("ERROR:  OUTER RADIUS must be greater than 1.\n");
    exit(-1);
  }

  options.getArgs("RADIAL EXPANSION MODES", shell->Nmodes);

  if (options.compareArgs("RADIAL BASIS TYPE", "TRUE")) {
    shellLoadRadialBasisTrue(shell);
  } else if (options.compareArgs("RADIAL BASIS TYPE", "GLOBALDISCRETE")) {
    shellLoadRadialBasisGlobalDiscrete(shell);
  } else if (options.compareArgs("RADIAL BASIS TYPE", "PIECEWISEDISCRETE")) {
    shellLoadRadialBasisPiecewiseDiscrete(shell);
  } else {
    printf("ERROR:  Unrecognized value for option RADIAL BASIS TYPE.\n");
    exit(-1);
  }

  // Put the r^2 factor into the quadrature weights.
  //
  // TODO:  Fix other codes so this can be removed.
  for (int i = 0; i < shell->Nquad; i++)
    shell->Wquad[i] *= pow(shell->Rquad[i], 2);

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
  shell->Ntotal = Nlocal + Nhalo;

  shell->r3D = (dfloat*) calloc(shell->Nmodes*shell->Ntotal, sizeof(dfloat));
  shell->q3D = (dfloat*) calloc(shell->Nmodes*shell->Ntotal, sizeof(dfloat));
  shell->f   = (dfloat*) calloc(shell->Nquad, sizeof(dfloat));
  shell->r   = (dfloat*) calloc(shell->Ntotal, sizeof(dfloat));
  shell->x   = (dfloat*) calloc(shell->Ntotal, sizeof(dfloat));

  kernelInfo["defines/p_shellNmodes"] = shell->Nmodes;
  kernelInfo["defines/p_shellNquad"] = shell->Nquad;
  kernelInfo["defines/p_shellNgll"] = shell->Ngll;

  occa::properties kernelInfoP  = kernelInfo;

  shell->o_r   = mesh->device.malloc(shell->Ntotal*sizeof(dfloat), shell->r);
  shell->o_x   = mesh->device.malloc(shell->Ntotal*sizeof(dfloat), shell->x);

  //Solver tolerances
  shell->pTOL = 1E-8;

  shell->elliptic = new elliptic_t;
  shell->elliptic->mesh = mesh;
  shell->elliptic->options = shell->options;
  shell->elliptic->dim = shell->dim;
  shell->elliptic->elementType = shell->elementType;

  // Elliptic BCType flags should be all zero for spherical problem.
  shell->elliptic->BCType = (int*)calloc(7, sizeof(int));

  //ellipticSolveSetup(shell->elliptic, shell->lambda, kernelInfoP);

  shell->meshSEM = NULL;

  // TODO:  Make this work for Nradelements > 1.
  if ((shell->elementType == QUADRILATERALS) &&
      (options.compareArgs("RADIAL BASIS TYPE", "PIECEWISEDISCRETE")) &&
      (shell->Nradelements == 1)) {
    shellExtrudeSphere(shell);
  }

  // OKL kernels specific to shell
  shell->shellReconstructKernel =
    mesh->device.buildKernel(DSHELL "/okl/shellReconstructHex3D.okl",
        "shellReconstructHex3D",
        kernelInfo);
}

/*****************************************************************************/

static void shellLoadRadialBasisTrue(shell_t *shell)
{
  dfloat *basisR, *basisLambda;    // Values of R, lambda for the cached basis.
  int Nrows, Ncols;                // Needed by readDfloatArray().
  FILE *fp;                        // Pointer to open data file.
  char fname[BUFSIZ];              // Path to data file.

  sprintf(fname, DHOLMES "/solvers/shell/data/shellN%02d.dat", shell->Nmodes);
  fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  readDfloatArray(fp, "SHELL TRUE R",
      &basisR, &Nrows, &Ncols);
  readDfloatArray(fp, "SHELL TRUE LAMBDA",
      &basisLambda, &Ncols, &Ncols);

  // Make sure the parameters match those corresponding to the cached basis.
  //
  // TODO:  These checks might need tolerances to avoid floating-point issues.
  if ((shell->lambda != *basisLambda) || (shell->R != *basisR)) {
    printf("ERROR:  Basis R = %f, lambda = %f; got R = %f, lambda = %f.\n",
           *basisR, *basisLambda, shell->R, shell->lambda);
    printf("Please recompute the true eigenfunction basis.");
    exit(-1);
  }

  // Load the basis data.
  if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
    // Dirichlet BCs on both spheres.
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE NODES",
        &(shell->Rquad), &(shell->Nquad), &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE WEIGHTS",
        &(shell->Wquad), &(shell->Nquad), &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET GLL NODES",
        &(shell->Rgll), &(shell->Ngll), &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET PLOT NODES",
        &(shell->Rplot), &(shell->Nplot), &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET EIGENVALUES",
        &(shell->eigenvalues), &(shell->Nmodes), &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE",
        &(shell->Bquad), &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE",
        &(shell->DBquad), &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET GLL VANDERMONDE",
        &(shell->Bgll), &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE",
        &(shell->DBgll), &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET PLOT VANDERMONDE",
        &(shell->Bplot), &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL TRUE DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE",
        &(shell->DBplot), &Nrows, &Ncols);
  } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
    // Neumann BCs on both spheres.
    printf("Neumann-Neumann boundary conditions not implemented.\n");
    exit(-1);
  } else {
    printf("Type-{%d, %d} boundary conditions not implmented.\n", 
           shell->innerBC, shell->outerBC);
    exit(-1);
  }

  // We're done. Clean up.
  free(basisR);
  free(basisLambda);

  fclose(fp);
}

static void shellLoadRadialBasisGlobalDiscrete(shell_t *shell)
{
  dfloat *Tquad, *Tgll, *Tplot;    // Initial basis Vandermonde matrices.
  dfloat *DTquad, *DTgll, *DTplot; // Initial basis derivative Vandermondes.
  dfloat *A, *B;                   // Matrices for GEVP.
  dfloat C1, C2;                   // Scaling constants (e.g., Jacobian).
  int Nrows, Ncols;                // Needed by readDfloatArray().
  int Nmodes, Nquad, Ngll, Nplot;  // Local copies of variables from shell.
  FILE *fp;                        // Pointer to open data file.
  char fname[BUFSIZ];              // Path to data file.

  // Variables needed for dsygv_()---see LAPACK documentation.
  int ITYPE, N, LDA, LDB, LWORK, INFO;
  char JOBZ, UPLO;
  dfloat *W, *WORK;

  sprintf(fname, DHOLMES "/solvers/shell/data/shellN%02d.dat", shell->Nmodes);
  fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  // Load the various node sets and quadrature weights.
  readDfloatArray(fp, "SHELL GLOBAL DISCRETE QUADRATURE NODES",
      &(shell->Rquad),&(shell->Nquad), &Ncols);
  readDfloatArray(fp, "SHELL GLOBAL DISCRETE QUADRATURE WEIGHTS",
      &(shell->Wquad),&(shell->Nquad), &Ncols);
  readDfloatArray(fp, "SHELL GLOBAL DISCRETE GLL NODES",
      &(shell->Rgll),&(shell->Ngll), &Ncols);
  readDfloatArray(fp, "SHELL GLOBAL DISCRETE PLOT NODES",
      &(shell->Rplot),&(shell->Nplot), &Ncols);

  // Scale nodes and weights from [-1, 1] to [1 R].
  for (int i = 0; i < shell->Nquad; i++) {
    shell->Rquad[i] = (1 + shell->Rquad[i])*((shell->R - 1)/2) + 1;
    shell->Wquad[i] = shell->Wquad[i]*(shell->R - 1)/2;
  }

  for (int i = 0; i < shell->Ngll; i++)
    shell->Rgll[i] = (1 + shell->Rgll[i])*((shell->R - 1)/2) + 1;
  for (int i = 0; i < shell->Nplot; i++)
    shell->Rplot[i] = (1 + shell->Rplot[i])*((shell->R - 1)/2) + 1;

  Nmodes = shell->Nmodes;
  Nquad  = shell->Nquad;
  Ngll   = shell->Ngll;
  Nplot  = shell->Nplot;

  // Load the initial basis matrices.
  if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
    // Dirichlet BCs on both spheres.
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE",
        &Tquad, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE",
        &DTquad, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL VANDERMONDE",
        &Tgll, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE",
        &DTgll, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT VANDERMONDE",
        &Tplot, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE",
        &DTplot, &Nrows, &Ncols);

    C1 = pow((shell->R - 1)/2, 2);
    C2 = (shell->R - 1)/2;
  } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
    // Neumann BCs on both spheres.
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE VANDERMONDE",
        &Tquad, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE DERIVATIVE VANDERMONDE",
        &DTquad, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL VANDERMONDE",
        &Tgll, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL DERIVATIVE VANDERMONDE",
        &DTgll, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT VANDERMONDE",
        &Tplot, &Nrows, &Ncols);
    readDfloatArray(fp, "SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT DERIVATIVE VANDERMONDE",
        &DTplot, &Nrows, &Ncols);

    C1 = pow((shell->R - 1), 3)/8;
    C2 = pow((shell->R - 1)/2, 2);
  } else {
    printf("Type-{%d, %d} boundary conditions not implmented.\n", 
           shell->innerBC, shell->outerBC);
    exit(-1);
  }

  // Scale initial basis matrices from [-1, 1] to [1, R].
  for (int i = 0; i < Nquad*Nmodes; i++) {
      Tquad[i] *= C1;
      DTquad[i] *= C2;
  }

  for (int i = 0; i < Ngll*Nmodes; i++) {
      Tgll[i] *= C1;
      DTgll[i] *= C2;
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
        a1 = shell->Wquad[k]*pow(shell->Rquad[k], 2)*DTquad[i + k*Nmodes]*DTquad[j + k*Nmodes];
        a2 = shell->lambda*shell->Wquad[k]*pow(shell->Rquad[k], 2)*Tquad[i + k*Nmodes]*Tquad[j + k*Nmodes];
        b = shell->Wquad[k]*Tquad[i + k*Nmodes]*Tquad[j + k*Nmodes];
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
  shell->eigenvalues = (dfloat*) calloc(Nmodes, sizeof(dfloat));
  memcpy(shell->eigenvalues, W, Nmodes*sizeof(dfloat));

  // TODO:  Explicitly normalize the eigenfunctions?

  // Assemble Vandermonde matrices.
  //
  // NB:  A contains the eigenvectors, stored in column-major order.
  shell->Bquad = (dfloat*) calloc(Nquad*Nmodes, sizeof(dfloat));
  shell->DBquad = (dfloat*) calloc(Nquad*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Nquad; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        shell->Bquad[i*Nmodes + j]  += Tquad[i*Nmodes + k]*A[j*Nmodes + k];
        shell->DBquad[i*Nmodes + j] += DTquad[i*Nmodes + k]*A[j*Nmodes + k];
      }
    }
  }

  shell->Bgll = (dfloat*) calloc(Ngll*Nmodes, sizeof(dfloat));
  shell->DBgll = (dfloat*) calloc(Ngll*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Ngll; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        shell->Bgll[i*Nmodes + j]  += Tgll[i*Nmodes + k]*A[j*Nmodes + k];
        shell->DBgll[i*Nmodes + j] += DTgll[i*Nmodes + k]*A[j*Nmodes + k];
      }
    }
  }

  shell->Bplot = (dfloat*) calloc(Nplot*Nmodes, sizeof(dfloat));
  shell->DBplot = (dfloat*) calloc(Nplot*Nmodes, sizeof(dfloat));
  for (int i = 0; i < Nplot; i++) {
    for (int j = 0; j < Nmodes; j++) {
      for (int k = 0; k < Nmodes; k++) {
        shell->Bplot[i*Nmodes + j]   += Tplot[i*Nmodes + k]*A[j*Nmodes + k];
        shell->DBplot[i*Nmodes + j]  += DTplot[i*Nmodes + k]*A[j*Nmodes + k];
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

  fclose(fp);
}


// WARNING:  This will not do the right thing if Nradelements == 1.
//
// TODO:  Fix this.
static void shellLoadRadialBasisPiecewiseDiscrete(shell_t *shell)
{
  dfloat *Rquadb, *Wquadb;    // Base quadrature nodes and weights.
  dfloat *Rplotb;             // Base plotting nodes.
  dfloat *Rgllb;              // Base GLL nodes.
  dfloat *Bquadb, *DBquadb;   // Base Vandermonde matrices.
  dfloat *Bplotb, *DBplotb;
  dfloat r, J;                // Radial variable and Jacobian factor.
  dfloat *A, *B;              // Matrices for GEVP.
  int N, Nradelements, Nqr;   // Local copies of variables from shell.
  int Nplotr, Ngllr, Nmodes;
  int Nrows, Ncols;           // Needed by readDfloatArray().
  int ee, es, ind, off, end;  // Variables to assist with indexing.
  FILE *fp;                   // Pointer to open data file.
  char fname[BUFSIZ];         // Path to data file.

  // Variables needed for dsygv_()---see LAPACK documentation.
  int ITYPE, M, LDA, LDB, LWORK, INFO;
  char JOBZ, UPLO;
  dfloat *W, *WORK;

  N = shell->mesh->N;

  sprintf(fname, DHOLMES "/solvers/shell/data/shellN%02d.dat", N);
  fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  // Load base sets of quadrature nodes and weights and Vandermonde matrices.
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE QUADRATURE NODES",
      &Rquadb, &(shell->Nqr), &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE QUADRATURE WEIGHTS",
      &Wquadb, &(shell->Nqr), &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE PLOT NODES",
      &Rplotb, &(shell->Nplotr), &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE GLL NODES",
      &Rgllb, &(shell->Ngllr), &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE QUADRATURE VANDERMONDE",
      &Bquadb, &Nrows, &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE QUADRATURE DERIVATIVE VANDERMONDE",
      &DBquadb, &Nrows, &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE PLOT VANDERMONDE",
      &Bplotb, &Nrows, &Ncols);
  readDfloatArray(fp, "SHELL PIECEWISE DISCRETE PLOT DERIVATIVE VANDERMONDE",
      &DBplotb, &Nrows, &Ncols);

  Nradelements = shell->Nradelements;
  Nqr          = shell->Nqr;
  Nplotr       = shell->Nplotr;
  Ngllr        = shell->Ngllr;

  // Scale the base nodes and weights to get nodes and weights for each radial
  // subinterval.
  shell->Nquad = Nqr*Nradelements;
  shell->Nplot = Nplotr*Nradelements;
  shell->Ngll  = Ngllr*Nradelements;
  shell->Rquad = (dfloat*)calloc(shell->Nquad, sizeof(dfloat));
  shell->Wquad = (dfloat*)calloc(shell->Nquad, sizeof(dfloat));
  shell->Rplot = (dfloat*)calloc(shell->Nplot, sizeof(dfloat));
  shell->Rgll  = (dfloat*)calloc(shell->Ngll, sizeof(dfloat));
  for (int e = 0; e < Nradelements; e++) {
    J = (shell->Rbreaks[e + 1] - shell->Rbreaks[e])/2.0;

    for (int i = 0; i < Nqr; i++) {
      shell->Rquad[e*Nqr + i] = (Rquadb[i] + 1.0)*J + shell->Rbreaks[e];
      shell->Wquad[e*Nqr + i] = Wquadb[i]*J;
    }

    for (int i = 0; i < Nplotr; i++)
      shell->Rplot[e*Nplotr + i] = (Rplotb[i] + 1.0)*J + shell->Rbreaks[e];

    for (int i = 0; i < Ngllr; i++)
      shell->Rgll[e*Ngllr + i] = (Rgllb[i] + 1.0)*J + shell->Rbreaks[e];
  }

  // Determine the number of modes and set up some indexing variables.
  if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
    // Dirichlet-Dirichlet
    shell->Nmodes = Nradelements*N - 1;
  } else if (((shell->innerBC == 1) && (shell->outerBC == 2)) ||
             ((shell->innerBC == 2) && (shell->outerBC == 1))) {
    // Dirichlet-Neumann or Neumann-Dirichlet
    shell->Nmodes = Nradelements*N;
  } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
    // Neumann-Neumann
    shell->Nmodes = Nradelements*N + 1;
  }

  Nmodes = shell->Nmodes;

  if (shell->innerBC == 1) {
    es = 1;
  } else {
    es = 0;
  }

  if (shell->outerBC == 1) {
    ee = Nradelements - 1;
  } else {
    ee = Nradelements;
  }

  // Set up the eigenvalue problem.
  A = (dfloat*)calloc(Nmodes*Nmodes, sizeof(dfloat));
  B = (dfloat*)calloc(Nmodes*Nmodes, sizeof(dfloat));

  if (Nradelements == 1) {
    if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
      off = 1;
      end = N - 1;
    } else if ((shell->innerBC == 1) && (shell->outerBC == 2)) {
      off = 1;
      end = N;
    } else if ((shell->innerBC == 2) && (shell->outerBC == 1)) {
      off = 0;
      end = N;
    } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
      off = 0;
      end = N + 1;
    }

    J = (shell->Rbreaks[1] - shell->Rbreaks[0])/2.0;
    for (int i = 0; i < end; i++) {
      for (int j = 0; j < end; j++) {
        for (int k = 0; k < Nqr; k++) {
          ind = i*Nmodes + j;
          r = (Rquadb[k] + 1.0)*J + shell->Rbreaks[0];
          A[ind] += r*r*DBquadb[k*(N + 1) + i + off]*DBquadb[k*(N + 1) + j + off]*Wquadb[k]/J;
          A[ind] += shell->lambda*r*r*Bquadb[k*(N + 1) + i + off]*Bquadb[k*(N + 1) + j + off]*Wquadb[k]*J;
          B[ind] += Bquadb[k*(N + 1) + i + off]*Bquadb[k*(N + 1) + j + off]*Wquadb[k]*J;
        }
      }
    }
  } else {
    if (shell->innerBC == 1) {
      off = 0;
      J = (shell->Rbreaks[1] - shell->Rbreaks[0])/2.0;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < Nqr; k++) {
            ind = (i + off)*Nmodes + (j + off);
            r = (Rquadb[k] + 1.0)*J + shell->Rbreaks[0];
            A[ind] += r*r*DBquadb[k*(N + 1) + i + 1]*DBquadb[k*(N + 1) + j + 1]*Wquadb[k]/J;
            A[ind] += shell->lambda*r*r*Bquadb[k*(N + 1) + i + 1]*Bquadb[k*(N + 1) + j + 1]*Wquadb[k]*J;
            B[ind] += Bquadb[k*(N + 1) + i + 1]*Bquadb[k*(N + 1) + j + 1]*Wquadb[k]*J;
          }
        }
      }
    }

    for (int e = es; e < ee; e++) {
      if (shell->innerBC == 1) {
        off = e*N - 1;
      } else {
        off = e*N;
      }

      J = (shell->Rbreaks[e + 1] - shell->Rbreaks[e])/2.0;
      for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N + 1; j++) {
          for (int k = 0; k < Nqr; k++) {
            ind = (i + off)*Nmodes + (j + off);
            r = (Rquadb[k] + 1.0)*J + shell->Rbreaks[e];
            A[ind] += r*r*DBquadb[k*(N + 1) + i]*DBquadb[k*(N + 1) + j]*Wquadb[k]/J;
            A[ind] += shell->lambda*r*r*Bquadb[k*(N + 1) + i]*Bquadb[k*(N + 1) + j]*Wquadb[k]*J;
            B[ind] += Bquadb[k*(N + 1) + i]*Bquadb[k*(N + 1) + j]*Wquadb[k]*J;
          }
        }
      }
    }

    if (shell->outerBC == 1) {
      off = Nmodes - N;
      J = (shell->Rbreaks[Nradelements - 1] - shell->Rbreaks[Nradelements - 2])/2.0;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          for (int k = 0; k < Nqr; k++) {
            ind = (i + off)*Nmodes + (j + off);
            r = (Rquadb[k] + 1.0)*J + shell->Rbreaks[Nradelements - 1];
            A[ind] += r*r*DBquadb[k*(N + 1) + i]*DBquadb[k*(N + 1) + j]*Wquadb[k]/J;
            A[ind] += shell->lambda*r*r*Bquadb[k*(N + 1) + i]*Bquadb[k*(N + 1) + j]*Wquadb[k]*J;
            B[ind] += Bquadb[k*(N + 1) + i]*Bquadb[k*(N + 1) + j]*Wquadb[k]*J;
          }
        }
      }
    }
  }

  // Solve the eigenvalue problem.
  ITYPE = 1;        // Solve Ax = cBx
  JOBZ = 'V';       // Want both eigenvalues and eigenvectors.
  UPLO = 'U';       // Assume upper triangular storage for symmetric matrices.
  M = Nmodes;       // Matrix dimension.
  LDA = Nmodes;     // Leading dimension of A.
  LDB = Nmodes;     // Leading dimension of B.
  LWORK = 4*Nmodes; // Workspace size.  (TODO:  ilaenv() for optimal value?)

  W = (dfloat*)calloc(M, sizeof(dfloat));        // Eigenvalues.
  WORK = (dfloat*)calloc(LWORK, sizeof(dfloat)); // LAPACK workspace.

  dsygv_(&ITYPE, &JOBZ, &UPLO, &M, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);

  if (INFO != 0)
    printf("Error in DSYGV (INFO = %d).\n", INFO);

  // Assign eigenvalues.  (TODO:  Sort them and the modes first?)
  shell->eigenvalues = (dfloat*)calloc(Nmodes, sizeof(dfloat));
  memcpy(shell->eigenvalues, W, Nmodes*sizeof(dfloat));

  // Assemble the Vandermonde matrices for the eigenmodes (global).
  //
  // NB:  A contains the eigenvectors, stored in column-major order.
  shell->Bquad = (dfloat*)calloc(shell->Nquad*Nmodes, sizeof(dfloat));
  shell->DBquad = (dfloat*)calloc(shell->Nquad*Nmodes, sizeof(dfloat));
  shell->Bplot = (dfloat*)calloc(shell->Nplot*Nmodes, sizeof(dfloat));
  shell->DBplot = (dfloat*)calloc(shell->Nplot*Nmodes, sizeof(dfloat));

  if (Nradelements == 1) {
    if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
      off = 1;
      end = N - 1;
    } else if ((shell->innerBC == 1) && (shell->outerBC == 2)) {
      off = 1;
      end = N;
    } else if ((shell->innerBC == 2) && (shell->outerBC == 1)) {
      off = 0;
      end = N;
    } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
      off = 0;
      end = N + 1;
    }

    J = (shell->Rbreaks[1] - shell->Rbreaks[0])/2.0;
    for (int i = 0; i < shell->Nqr; i++) {
      for (int j = 0; j < Nmodes; j++) {
        for (int k = 0; k < end; k++) {
          ind = i*shell->Nmodes + j;
          shell->Bquad[ind] += Bquadb[i*(N + 1) + k + off]*A[j*shell->Nmodes + k];
          shell->DBquad[ind] += DBquadb[i*(N + 1) + k + off]*A[j*shell->Nmodes + k]/J;
        }
      }
    }

    for (int i = 0; i < shell->Nplotr; i++) {
      for (int j = 0; j < Nmodes; j++) {
        for (int k = 0; k < end; k++) {
          ind = i*shell->Nmodes + j;
          shell->Bplot[ind] += Bplotb[i*(N + 1) + k + off]*A[j*shell->Nmodes + k];
          shell->DBplot[ind] += DBplotb[i*(N + 1) + k + off]*A[j*shell->Nmodes + k]/J;
        }
      }
    }
  } else {
    if (shell->innerBC == 1) {
      J = (shell->Rbreaks[1] - shell->Rbreaks[0])/2.0;
      for (int i = 0; i < shell->Nqr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N; k++) {
            ind = i*shell->Nmodes + j;
            shell->Bquad[ind] += Bquadb[i*(N + 1) + k + 1]*A[j*shell->Nmodes + k];
            shell->DBquad[ind] += DBquadb[i*(N + 1) + k + 1]*A[j*shell->Nmodes + k]/J;
          }
        }
      }

      for (int i = 0; i < shell->Nplotr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N; k++) {
            ind = i*shell->Nmodes + j;
            shell->Bplot[ind] += Bplotb[i*(N + 1) + k + 1]*A[j*shell->Nmodes + k];
            shell->DBplot[ind] += DBplotb[i*(N + 1) + k + 1]*A[j*shell->Nmodes + k]/J;
          }
        }
      }

      off = -1;
    } else {
      off = 0;
    }

    for (int e = es; e < ee; e++) {
      J = (shell->Rbreaks[e + 1] - shell->Rbreaks[e])/2.0;

      for (int i = 0; i < shell->Nqr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N + 1; k++) {
            ind = (i + e*shell->Nqr)*shell->Nmodes + j;
            shell->Bquad[ind] += Bquadb[i*(N + 1) + k]*A[j*shell->Nmodes + e*N + k + off];
            shell->DBquad[ind] += DBquadb[i*(N + 1) + k]*A[j*shell->Nmodes + e*N + k + off]/J;
          }
        }
      }

      for (int i = 0; i < shell->Nplotr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N + 1; k++) {
            ind = (i + e*shell->Nplotr)*shell->Nmodes + j;
            shell->Bplot[ind] += Bplotb[i*(N + 1) + k]*A[j*shell->Nmodes + e*N + k + off];
            shell->DBplot[ind] += DBplotb[i*(N + 1) + k]*A[j*shell->Nmodes + e*N + k + off]/J;
          }
        }
      }
    }

    if (shell->outerBC == 1) {
      J = (shell->Rbreaks[Nradelements - 1] - shell->Rbreaks[Nradelements - 2])/2.0;

      for (int i = 0; i < shell->Nqr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N; k++) {
            ind = (i + (Nradelements - 1)*shell->Nqr)*shell->Nmodes + j;
            shell->Bquad[ind] += Bquadb[i*(N + 1) + k]*A[j*shell->Nmodes + (Nradelements - 1)*N + k + off];
            shell->DBquad[ind] += DBquadb[i*(N + 1) + k]*A[j*shell->Nmodes + (Nradelements - 1)*N + k + off]/J;
          }
        }
      }

      for (int i = 0; i < shell->Nplotr; i++) {
        for (int j = 0; j < Nmodes; j++) {
          for (int k = 0; k < N; k++) {
            ind = (i + (Nradelements - 1)*shell->Nplotr)*shell->Nmodes + j;
            shell->Bplot[ind] += Bplotb[i*(N + 1) + k]*A[j*shell->Nmodes + (Nradelements - 1)*N + k + off];
            shell->DBplot[ind] += DBplotb[i*(N + 1) + k]*A[j*shell->Nmodes + (Nradelements - 1)*N + k + off]/J;
          }
        }
      }
    }
  }

  // We're done.  Clean up.
  free(A);
  free(B);
  free(W);
  free(WORK);
  free(Rquadb);
  free(Wquadb);
  free(Rplotb);
  free(Rgllb);
  free(Bquadb);
  free(DBquadb);
  free(Bplotb);
  free(DBplotb);

  fclose(fp);
}
