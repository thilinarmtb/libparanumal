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

static void stokesTestInitialGuessZero2D(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p);

static void stokesTestSolutionDebug(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p);
static void stokesTestInitialGuessDebug(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p);
static void stokesTestForcingFunctionDebug(dfloat x, dfloat y, dfloat lambda, dfloat *fx, dfloat *fy);

static void stokesTestSolutionSimpleConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p);
static void stokesTestForcingFunctionSimpleConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *fx, dfloat *fy);

static void stokesTestSolutionSimpleTimeDependentConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p);
static void stokesTestForcingFunctionSimpleTimeDependentConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *fx, dfloat *fy);

static void stokesTestSolutionSimpleConstantViscosityDirichletHex3D(dfloat x, dfloat y, dfloat z, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *uz, dfloat *p);
static void stokesTestForcingFunctionSimpleConstantViscosityDirichletHex3D(dfloat x, dfloat y, dfloat z, dfloat lambda, dfloat *fx, dfloat *fy, dfloat *fz);

/* TODO: Accommodate pressure forcing?  (Changes to the RHS code needed.) */

/* TODO:  Support variable viscosity. */

void stokesGetTestCase(stokes_t *stokes, stokesTestCase_t *testCase)
{
  string name;

  stokes->options.getArgs("TEST CASE", name);

  testCase->isTimeDependent = 0;
  testCase->solFn2D         = NULL;
  testCase->solFn3D         = NULL;
  testCase->initGuessFn2D   = NULL;
  testCase->initGuessFn3D   = NULL;
  testCase->forcingFn2D     = NULL;
  testCase->forcingFn3D     = NULL;
  testCase->tdSolFn2D       = NULL;
  testCase->tdSolFn3D       = NULL;
  testCase->tdForcingFn2D   = NULL;
  testCase->tdForcingFn3D   = NULL;

  /* TODO:  Check the dimension here to save ourselves from some stupid mistakes. */
  if (stokes->meshV->dim == 2) {
    if (name == "Debug") {
      testCase->solFn2D       = stokesTestSolutionDebug;
      testCase->initGuessFn2D = stokesTestInitialGuessDebug;
      //testCase->initGuessFn2D = stokesTestInitialGuessZero2D;
      testCase->forcingFn2D   = stokesTestForcingFunctionDebug;
    } else if (name == "SimpleConstantViscosityDirichletQuad2D") {
      testCase->solFn2D     = stokesTestSolutionSimpleConstantViscosityDirichletQuad2D;
      testCase->forcingFn2D = stokesTestForcingFunctionSimpleConstantViscosityDirichletQuad2D;
    } else if (name == "SimpleTimeDependentConstantViscosityDirichletQuad2D") {
      testCase->isTimeDependent = 1;
      testCase->tdSolFn2D     = stokesTestSolutionSimpleTimeDependentConstantViscosityDirichletQuad2D;
      testCase->tdForcingFn2D = stokesTestForcingFunctionSimpleTimeDependentConstantViscosityDirichletQuad2D;
    } else {
      printf("ERROR:  Invalid 2D test case %s.\n", name.c_str());
      MPI_Finalize();
      exit(-1);
    }
  } else if (stokes->meshV->dim == 3) {
    if (name == "SimpleConstantViscosityDirichletHex3D") {
      testCase->solFn3D     = stokesTestSolutionSimpleConstantViscosityDirichletHex3D;
      testCase->forcingFn3D = stokesTestForcingFunctionSimpleConstantViscosityDirichletHex3D;
    } else {
      printf("ERROR:  Invalid 3D test case %s.\n", name.c_str());
      MPI_Finalize();
      exit(-1);
    }
  }

  return;
}

static void stokesTestInitialGuessZero2D(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p)
{
  *ux = 0.0;
  *uy = 0.0;
  *p = 0.0;
  return;
}

/* TODO:  Need to setup kernels for boundary conditions. */

/*****************************************************************************/
/* Debug
 *
 * Scratch test case for debugging.
 */

static void stokesTestSolutionDebug(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p)
{
  const dfloat t = 1.0e-1;
  *ux = cos(y)*cos(t);
  *uy = sin(x)*sin(t);
  *p  = x*cos(t) + y*sin(t);
  return;
}

static void stokesTestInitialGuessDebug(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p)
{
  *ux = cos(y);
  *uy = 0.0;
  *p = x;
  return;
}

static void stokesTestForcingFunctionDebug(dfloat x, dfloat y, dfloat lambda, dfloat *fx, dfloat *fy)
{
  const dfloat t = 1.0e-1;
  *fx = cos(y)/1.0e-1 + cos(t) + (cos(t) - sin(t))*cos(y);
  *fy = 0.0           + sin(t) + (cos(t) + sin(t))*sin(x);
  return;
}

/*****************************************************************************/
/* SimpleConstantViscosityDirichletQuad2D
 *
 * Exact solution:  u_x = cos(y)
 *                  u_y = sin(x)
 *                  p   = x + y
 */

static void stokesTestSolutionSimpleConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *p)
{
  *ux = cos(y);
  *uy = sin(x);
  *p  = x + y;
  return;
}

static void stokesTestForcingFunctionSimpleConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat lambda, dfloat *fx, dfloat *fy)
{
  *fx = 1.0 + (1.0 + lambda)*cos(y);
  *fy = 1.0 + (1.0 + lambda)*sin(x);
  return;
}

/*****************************************************************************/
/* SimpleTimeDependentConstantViscosityDirichletQuad2D
 *
 * Exact solution:  u_x = cos(y)cos(t)
 *                  u_y = sin(x)sin(t)
 *                  p   = x*cos(t) + y*sin(t)
 */

static void stokesTestSolutionSimpleTimeDependentConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat t, dfloat *ux, dfloat *uy, dfloat *p)
{
  *ux = cos(y)*cos(t);
  *uy = sin(x)*sin(t);
  *p  = x*cos(t) + y*sin(t);
  return;
}

static void stokesTestForcingFunctionSimpleTimeDependentConstantViscosityDirichletQuad2D(dfloat x, dfloat y, dfloat t, dfloat *fx, dfloat *fy)
{
  *fx = cos(t) + (cos(t) - sin(t))*cos(y);
  *fy = sin(t) + (cos(t) + sin(t))*sin(x);
  return;
}

/*****************************************************************************/
/* SimpleConstantViscosityDirichletHex3D
 *
 * Exact solution:  u_x = -6z(1 - z^2)^2
 *                  u_y = -6x(1 - x^2)^2
 *                  u_z = -6y(1 - y^2)^2
 *                  p   = sin(pi*x)sin(pi*y)sin(pi*z)
 */

static void stokesTestSolutionSimpleConstantViscosityDirichletHex3D(dfloat x, dfloat y, dfloat z, dfloat lambda, dfloat *ux, dfloat *uy, dfloat *uz, dfloat *p)
{
  *ux = -6.0*z*pow(1.0 - z*z, 2.0);
  *uy = -6.0*x*pow(1.0 - x*x, 2.0);
  *uz = -6.0*y*pow(1.0 - y*y, 2.0);
  *p  = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
  return;
}

static void stokesTestForcingFunctionSimpleConstantViscosityDirichletHex3D(dfloat x, dfloat y, dfloat z, dfloat lambda, dfloat *fx, dfloat *fy, dfloat *fz)
{
  *fx = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 48.0*z*z*z - 72.0*z*(1.0 - z*z) - lambda*6.0*z*pow(1.0 - z*z, 2.0);
  *fy = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z) + 48.0*x*x*x - 72.0*x*(1.0 - x*x) - lambda*6.0*x*pow(1.0 - x*x, 2.0);
  *fz = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z) + 48.0*y*y*y - 72.0*y*(1.0 - y*y) - lambda*6.0*y*pow(1.0 - y*y, 2.0);
  return;
}

/* TODO:  Reactivate these test cases. */

#if 0

/* All-Neumann BC test case. */
static void stokesTestSolutionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy)
{
  *ux = 6.0*pow(1.0 - x*x, 3.0)*pow(1.0 - y*y, 2.0)*y;
  *uy = -6.0*pow(1.0 - y*y, 3.0)*pow(1.0 - x*x, 2.0)*x;
  return;
}

static void stokesTestForcingFunctionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = cos(M_PI*x)*sin(M_PI*y)/M_PI - 24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) - 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0);
  *fy = sin(M_PI*x)*cos(M_PI*y)/M_PI + 24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0);
  return;
}

/* Variable viscosity Quad2D test case with eta = 2.0 + sinh(x*y) */
static void stokesTestSolutionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy)
{
  *ux = 6.0*pow(1.0 - x*x, 3.0)*pow(1.0 - y*y, 2.0)*y;
  *uy = -6.0*pow(1.0 - y*y, 3.0)*pow(1.0 - x*x, 2.0)*x;
  return;
}

static void stokesTestForcingFunctionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = cos(M_PI*x)*sin(M_PI*y)/M_PI - (2.0 + sinh(x*y))*(24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) + 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0)) + 36.0*x*pow(1.0 - x*x, 2.0)*pow(1.0 - y*y, 2.0)*y*y*cosh(x*y) - 6.0*pow(1.0 - x*x, 3.0)*(1.0 - 5.0*y*y)*(1.0 - y*y)*x*cosh(x*y);
  *fy = sin(M_PI*x)*cos(M_PI*y)/M_PI + (2.0 + sinh(x*y))*(24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0)) + 6.0*pow(1.0 - y*y, 3.0)*(1.0 - 5.0*x*x)*(1.0 - x*x)*y*cosh(x*y) - 36.0*y*pow(1.0 - y*y, 2.0)*pow(1.0 - x*x, 2.0)*x*x*cosh(x*y);
  return;
}

/* Quad2D leaky cavity test case.  Need to set up BCs.  This one has no closed-form solution. */
static void stokesTestForcingFunctionLeakyCavityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = 0.0;
  *fy = 0.0;
  return;
}

#endif

