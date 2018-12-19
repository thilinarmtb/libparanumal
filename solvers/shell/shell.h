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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"

// For the ASBF shell solver, always use JACOBI for modes indexed this and
// higher.
//
// TODO:  This should be handled in the setup file.
#define SHELL_ASBF_JACOBI_CROSSOVER 6

typedef struct {

  int dim, elementType;

  mesh_t     *mesh;
  mesh_t     *meshSEM;  // 3D hex SEM mesh
  elliptic_t *elliptic;
  precon_t   *preconJacobi;
  
  setupAide options;

  /*********/
  dlong Ntotal;

  int Nmodes; // number of radial modes
  int Nquad;  // number of radial quadrature nodes
  int Ngll;   // number of radial gll nodes
  int Nplot;  // number of radial plotting points

  dfloat R;      // outer radius of shell; must have R > 1.

  dfloat *eigenvalues; // generalzied eigenvalues of discrete r^2 weighted 1D Laplacian on [1,1.5]

  dfloat *Rquad; // coordinates of radial quadrature nodes (start at 1)
  dfloat *Rgll;  // coordinates of radial gll nodes (start at 1)
  dfloat *Rplot; // coordinates of radial plot nodes (start at 1)
  
  dfloat *Wquad; // weights of radial quadrature nodes (include radius^2 factor)

  dfloat *Bquad; // Vandermonde for radial modes evaluated at quadrature nodes
  dfloat *Bgll;  // Vandermonde for radial modes evaluated at gll nodes
  dfloat *Bplot; // Vandermonde for radial modes evaluated at plot nodes

  dfloat *DBquad; // Vandermonde for derivative of radial modes evaluated at quadrature nodes
  dfloat *DBgll;  // Vandermonde for derivative of radial modes evaluated at gll nodes
  dfloat *DBplot; // Vandermonde for derivative of radial modes evaluated at plot nodes
  
  dfloat *r3D;    // right-hand sides for screened Poisson equation for each mode
  dfloat *q3D;    // solution to screened Poisson equation for each mode
  dfloat *f;      // forcing function

  // Boundary conditions on the inner and outer sphere.  Set to 1 or 2 for
  // Dirichlet or Neumann, respectively.
  int innerBC;
  int outerBC;

  // Variables for piecewise discrete basis.
  int Nradelements; // number of radial elements
  int Nqr;          // number of quadrature nodes on each radial element.
  int Nplotr;       // number of plot nodes on each radial element.
  int Ngllr;        // number of GLL nodes on each radial element.
  dfloat *Rbreaks;  // radial mesh breakpoints

  // Profiling information.
  struct {
    struct {
      double solveSetup;
      double rhsSetup;
      double total;
    } setup;

    double solve;

    double total;
  } times;

  /*********/

  dfloat lambda;    // helmhotz solver -lap(u) + lamda u

  //solver tolerances
  dfloat TOL;

  dfloat *r;
  dfloat *x;

  occa::memory o_r;
  occa::memory o_x;

  occa::kernel shellReconstructKernel;
  
} shell_t;

shell_t *shellSetup(mesh_t *mesh, dfloat lambda, occa::properties kernelInfo, setupAide options);

int shellSolve(shell_t *shell, setupAide options);

void shellSolveSetup(shell_t *shell, dfloat lambda, occa::properties &kernelInfo);

void shellErrorHex3D(shell_t *shell, dfloat *q, dfloat *normErrorL2, dfloat *normErrorH1);

void shellExtrudeSphere(shell_t *shell);

void shellPlotVTU3D(shell_t *shell, const char *fileNameBase, int fld);

void shellCubatureGradient(shell_t *shell, dfloat *q3D,
			  dfloat *cubq, dfloat *cubdqdx, dfloat *cubdqdy, dfloat *cubdqdz);
