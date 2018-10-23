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
  elliptic_t *pSolver;

  mesh_t *meshSEM; // 3D hex SEM mesh
  
  setupAide options;
  setupAide pOptions;

  /*********/
  dlong Ntotal;

  int asbfNmodes; // number of ASBF modes
  int asbfNquad;  // number of ASBF quadrature nodes
  int asbfNgll;   // number of ASBF gll nodes
  int asbfNplot;  // number of ASBF plotting points

  dfloat *asbfEigenvalues; // generalzied eigenvalues of discrete r^2 weighted 1D Laplacian on [1,1.5]

  dfloat *asbfRquad; // coordinates of ABSF quadrature nodes (start at 1)
  dfloat *asbfRgll;  // coordinates of ABSF gll nodes (start at 1)
  dfloat *asbfRplot; // coordinates of ABSF plot nodes (start at 1)
  
  dfloat *asbfWquad; // weights of ABSF quadrature nodes (include radius^2 factor)

  dfloat *asbfBquad; // generalized Vandermonde for ABSF modes evaluated at quadrature nodes
  dfloat *asbfBgll;  // generalized Vandermonde for ABSF modes evaluated at gll nodes
  dfloat *asbfBplot; // generalized Vandermonde for ABSF modes evaluated at plot nodes

  dfloat *asbfDBquad; // generalized Vandermonde for derivative of ABSF modes evaluated at quadrature nodes
  dfloat *asbfDBgll;  // generalized Vandermonde for derivative of ABSF modes evaluated at gll nodes
  dfloat *asbfDBplot; // generalized Vandermonde for derivative of ABSF modes evaluated at plot nodes
  
  dfloat *r3D;
  dfloat *q3D;
  dfloat *f;
  /*********/

  // ASBF SOLVER OCCA VARIABLES
  dfloat pbar;
  int Nblock;
  int frame;

  dfloat lambda;      // helmhotz solver -lap(u) + lamda u

  int NiterP;

  //solver tolerances
  dfloat pTOL;

  dfloat *r;
  dfloat *x;

  occa::memory o_r;
  occa::memory o_x;

  int *PmapB;
  occa::memory o_PmapB;

  //halo data
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;

  occa::memory o_pSendBuffer;
  occa::memory o_pRecvBuffer;

  occa::memory o_plotInterp, o_plotEToV;

  occa::kernel constrainKernel;

  occa::memory o_P;
  occa::memory o_rhsP;

  occa::memory o_pHaloBuffer;

  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;

  occa::kernel setFlowFieldKernel;

  occa::kernel pressureRhsKernel;
  occa::kernel pressureRhsIpdgBCKernel;
  occa::kernel pressureRhsBCKernel;
  occa::kernel pressureAddBCKernel;
  // occa::kernel pressureUpdateKernel;

  occa::kernel asbfReconstructKernel;
  
} asbf_t;

asbf_t *asbfSetup(mesh_t *mesh, setupAide options);

int asbfSolve(asbf_t *asbf, setupAide options);

void asbfErrorHex3D(asbf_t *asbf, dfloat *q);

void asbfExtrudeSphere(asbf_t *asbf);

void asbfPlotVTU3D(asbf_t *asbf, char *fileNameBase, int fld);

void asbfCubatureGradient(asbf_t *asbf, dfloat *q3D,
			  dfloat *cubq, dfloat *cubdqdx, dfloat *cubdqdy, dfloat *cubdqdz);

void interpolateQuad2D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);
