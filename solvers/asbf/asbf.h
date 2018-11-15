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

/* TODO:  THE ASBF DATA FIELDS NEED TO GO IN HERE. */
typedef struct {

  int dim, elementType;

  mesh_t *mesh;
  elliptic_t *elliptic;

  mesh_t *meshSEM; // 3D hex SEM mesh
  
  setupAide options;

  /*********/
  dlong Ntotal;

  int Nmodes; // number of ASBF modes
  int Nquad;  // number of ASBF quadrature nodes
  int Ngll;   // number of ASBF gll nodes
  int Nplot;  // number of ASBF plotting points

  dfloat R;      // outer radius of shell; must have R > 1.

  dfloat *eigenvalues; // generalzied eigenvalues of discrete r^2 weighted 1D Laplacian on [1,1.5]

  dfloat *Rquad; // coordinates of ABSF quadrature nodes (start at 1)
  dfloat *Rgll;  // coordinates of ABSF gll nodes (start at 1)
  dfloat *Rplot; // coordinates of ABSF plot nodes (start at 1)
  
  dfloat *Wquad; // weights of ABSF quadrature nodes (include radius^2 factor)

  dfloat *Bquad; // generalized Vandermonde for ABSF modes evaluated at quadrature nodes
  dfloat *Bgll;  // generalized Vandermonde for ABSF modes evaluated at gll nodes
  dfloat *Bplot; // generalized Vandermonde for ABSF modes evaluated at plot nodes

  dfloat *DBquad; // generalized Vandermonde for derivative of ABSF modes evaluated at quadrature nodes
  dfloat *DBgll;  // generalized Vandermonde for derivative of ABSF modes evaluated at gll nodes
  dfloat *DBplot; // generalized Vandermonde for derivative of ABSF modes evaluated at plot nodes
  
  dfloat *r3D;    // right-hand sides for screened Poisson equation for each mode
  dfloat *q3D;    // solution to screened Poisson equation for each mode
  dfloat *f;      // forcing function

  // Boundary conditions:  three integers {0, BCInner, BCOuter}.  Set BCInner
  // and BCOuter to 1 or 2 for Dirichlet or Neumann conditions, respectively.
  int *BCType;

  // Variables for piecewise discrete basis.
  int Nradelements; // number of radial elements
  int Nqr;          // number of quadrature nodes on each radial element.
  int Nplotr;       // number of plot nodes on each radial element.
  dfloat *Rbreaks;  // radial mesh breakpoints
  /*********/

  dfloat lambda;      // helmhotz solver -lap(u) + lamda u

  // ASBF SOLVER OCCA VARIABLES

  //solver tolerances
  dfloat pTOL;

  dfloat *r;
  dfloat *x;

  occa::memory o_r;
  occa::memory o_x;

  occa::kernel asbfReconstructKernel;
  
} asbf_t;

asbf_t *asbfSetup(mesh_t *mesh, dfloat lambda, occa::properties kernelInfo, setupAide options);

int asbfSolve(asbf_t *asbf, setupAide options);

void asbfSolveSetup(asbf_t *asbf, dfloat lambda, occa::properties &kernelInfo);

void asbfErrorHex3D(asbf_t *asbf, dfloat *q, dfloat *normErrorL2, dfloat *normErrorH1);

void asbfExtrudeSphere(asbf_t *asbf);

void asbfPlotVTU3D(asbf_t *asbf, char *fileNameBase, int fld);

void asbfCubatureGradient(asbf_t *asbf, dfloat *q3D,
			  dfloat *cubq, dfloat *cubdqdx, dfloat *cubdqdy, dfloat *cubdqdz);

void interpolateQuad2D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);
