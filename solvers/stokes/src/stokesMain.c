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

static void stokesTestSolutionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionDirichletQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionConstantViscosityHex3D(dfloat x, dfloat y, dfloat z, dfloat *ux, dfloat *uy, dfloat *uz);

int main(int argc, char **argv)
{
  dfloat           lambda;
  stokes_t         *stokes;
  stokesTestCase_t testCase;
  occa::properties kernelInfoV, kernelInfoP;

  // Start up MPI.
  MPI_Init(&argc, &argv);

  if (argc != 2) {
    printf("usage: ./stokesMain setupfile\n");
    MPI_Finalize();
    exit(-1);
  }

  setupAide options(argv[1]);

  options.getArgs("LAMBDA", lambda);

  stokes = stokesSetup(lambda, kernelInfoV, kernelInfoP, options);
  stokesSolve(stokes, lambda);

  stokesVecCopyDeviceToHost(stokes->u);

#if 1
  /* Compute error (if applicable.) */
  dfloat errxInf = 0.0, erryInf = 0.0, errzInf = 0.0;
  dfloat errxDL2 = 0.0, erryDL2 = 0.0, errzDL2 = 0.0;
  stokesGetTestCase(stokes, &testCase);
  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->Np; i++) {
      int    ind;
      dfloat x, y, z;
      dfloat errx, erry, errz;
      dfloat ux_exact, uy_exact, uz_exact, p_exact;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];
      z = stokes->meshV->z[ind];

      /* TODO:  Handle the case where the true solution is not known. */
      if (stokes->meshV->dim == 2) {
        stokesSolutionFunction2D solFn = (stokesSolutionFunction2D)testCase.solFn;
        solFn(x, y, lambda, &ux_exact, &uy_exact, &p_exact);

        /* Manually insert the boundary data.
         *
         * TODO:  Surely this should have been done already elsewhere?
         */
        if (stokes->mapB[e*stokes->meshV->Np + i] == 1) {
          stokes->u.x[ind] = ux_exact;
          stokes->u.y[ind] = uy_exact;
        }
      } else if (stokes->meshV->dim == 3) {
        stokesSolutionFunction3D solFn = (stokesSolutionFunction3D)testCase.solFn;
        solFn(x, y, z, lambda, &ux_exact, &uy_exact, &uz_exact, &p_exact);

        /* Manually insert the boundary data.
         *
         * TODO:  Surely this should have been done already elsewhere?
         */
        if (stokes->mapB[e*stokes->meshV->Np + i] == 1) {
          stokes->u.x[ind] = ux_exact;
          stokes->u.y[ind] = uy_exact;
          stokes->u.z[ind] = uz_exact;
        }
      }

      errx = stokes->u.x[ind] - ux_exact;
      erry = stokes->u.y[ind] - uy_exact;
      if (stokes->meshV->dim == 3)
        errz = stokes->u.z[ind] - uz_exact;

      if (fabs(errx) > errxInf)
        errxInf = fabs(errx);
      if (fabs(erry) > erryInf)
        erryInf = fabs(erry);
      errxDL2 += errx*errx;
      erryDL2 += erry*erry;

      if (stokes->meshV->dim == 3) {
        if (fabs(errz) > errzInf)
          errzInf = fabs(errz);
        errzDL2 += errz*errz;
      }
    }
  }

  errxDL2 = sqrt(errxDL2);
  erryDL2 = sqrt(erryDL2);
  if (stokes->meshV->dim == 3)
    errzDL2 = sqrt(errzDL2);

  printf("-----\n");

  printf("errxInf = % .15e\n", errxInf);
  printf("erryInf = % .15e\n", erryInf);
  if (stokes->meshV->dim == 3)
    printf("errzInf = % .15e\n", errzInf);
  printf("errxDL2 = % .15e\n", errxDL2);
  printf("erryDL2 = % .15e\n", erryDL2);
  if (stokes->meshV->dim == 3)
    printf("errzDL2 = % .15e\n", errzDL2);
#endif

  /* Export solution. */

  /* Report runtime statistics. */

  // Shut down MPI.
  MPI_Finalize();

  return 0;
}
