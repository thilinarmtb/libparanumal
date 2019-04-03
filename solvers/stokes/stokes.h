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

#ifndef STOKES_H
#define STOKES_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "parAlmond.hpp"
#include "elliptic.h"

/* Block size for reductions. */
#define STOKES_REDUCTION_BLOCK_SIZE 256

/* This data structure represents a vector used with the Stokes solver.  The
 * data is allocated in one big long block partitioned as
 *
 *   v = [ velocity x | velocity y | velocity z | pressure ].
 *
 * This is helpful, as we often don't care about the individual components but
 * about the vector as a whole.  For times when we want to address the
 * components individually, we also store pointers to the beginnings of each of
 * the component blocks for easy access.
 */
typedef struct {
  /* Host variables. */
  dfloat *v;          /* Vector data as one long block */
  dfloat *x;          /* Pointer to start of velocity x-component */
  dfloat *y;          /* Pointer to start of velocity y-component */
  dfloat *z;          /* Pointer to start of velocity z-component */
  dfloat *p;          /* Pointer to start of pressure */

  /* Device variables. */
  occa::memory o_v;
  occa::memory o_x;
  occa::memory o_y;
  occa::memory o_z;
  occa::memory o_p;
} stokesVec_t;

typedef struct {
  /* Jacobi preconditioner. */
  dfloat boost;          /* Jacobi boosting parameter. */
  stokesVec_t invDiagA;  /* Boosted inverse of the Stokes operator diagonal */


  /* Schur-complement-based block matrix preconditioner */
  elliptic_t *ellipticV; /* Elliptic solver object for velocity preconditioning */
  elliptic_t *ellipticP; /* Elliptic solver object for pressure multigrid */
  stokesVec_t invMM;     /* Boosted inverse of the pressure mass operator diagonal */
} stokesPrecon_t;

/* 2D test case / problem setup information. */
typedef void (*stokesSolutionFunction2D)(dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*);
typedef void (*stokesInitialGuessFunction2D)(dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*);
typedef void (*stokesForcingFunction2D)(dfloat, dfloat, dfloat, dfloat*, dfloat*);
typedef void (*stokesTimeDependentSolutionFunction2D)(dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*);
typedef void (*stokesTimeDependentForcingFunction2D)(dfloat, dfloat, dfloat, dfloat*, dfloat*);

/* 3D test case / problem setup information. */
typedef void (*stokesSolutionFunction3D)(dfloat, dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*, dfloat*);
typedef void (*stokesInitialGuessFunction3D)(dfloat, dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*, dfloat*);
typedef void (*stokesForcingFunction3D)(dfloat, dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*);
typedef void (*stokesTimeDependentSolutionFunction3D)(dfloat, dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*, dfloat*);
typedef void (*stokesTimeDependentForcingFunction3D)(dfloat, dfloat, dfloat, dfloat, dfloat*, dfloat*, dfloat*);

/* Data structure for holding information about test cases / problem setups. */
typedef struct {
  int                                   isTimeDependent;
  stokesSolutionFunction2D              solFn2D;
  stokesSolutionFunction3D              solFn3D;
  stokesInitialGuessFunction2D          initGuessFn2D;
  stokesInitialGuessFunction3D          initGuessFn3D;
  stokesForcingFunction2D               forcingFn2D;
  stokesForcingFunction3D               forcingFn3D;
  stokesTimeDependentSolutionFunction2D tdSolFn2D;
  stokesTimeDependentSolutionFunction3D tdSolFn3D;
  stokesTimeDependentForcingFunction2D  tdForcingFn2D;
  stokesTimeDependentForcingFunction3D  tdForcingFn3D;
} stokesTestCase_t;

typedef struct {
  setupAide options;       /* Configuration information. */

  int elementType;         /* Element type code */
  mesh_t *meshV;           /* Velocity mesh */
  mesh_t *meshP;           /* Pressure mesh */

  int NtotalV;             /* Total number of points in the velocity mesh */
  int NtotalP;             /* Total number of points in the pressure mesh */
  int Ndof;                /* Total number of degrees of freedom */

  stokesTestCase_t *testCase; /* Information/setup for the current test problem. */

  stokesVec_t u;           /* Solution */
  stokesVec_t f;           /* Right-hand side */

  dfloat      *eta;        /* Viscosity */
  occa::memory o_eta;

  dfloat      *cubEta;     /* Viscosity on cubature nodes */
  occa::memory o_cubEta;

  dfloat      *cubInterpV; /* Cubature interpolation matrix for the velocity mesh. */
  occa::memory o_cubInterpV;

  dfloat      *cubInterpP; /* Cubature interpolation matrix for the pressure mesh. */
  occa::memory o_cubInterpP;

  dfloat       *cubD;      /* Cubature differentiation matrix. */
  occa::memory o_cubD;

  stokesPrecon_t *precon;  /* Preconditioner */

  ogs_t *ogs;              /* Gather-scatter handle (masked for BCs) */

  /* Infrastructure for boundary conditions */
  int *BCType;             /* Physical-to-mathematical BC type code map. */
  dlong Nmasked;           /* Number of nodes masked out due to boundary conditions */

  dlong* maskIds;          /* Indices of the nodes to be masked */
  occa::memory o_maskIds;

  int* mapB;               /* Node-wise BC type codes */
  occa::memory o_mapB;

  /* infrastructure for Fischer successive RHS */
  int fsrNrhs;

  occa::memory o_fsrHistoryQ;
  occa::memory o_fsrZeroArray; // do not set entries of this
  occa::memory o_fsrAlphas;

  occa::kernel fsrStartKernel;
  occa::kernel fsrUpdateKernel;
  
  /* OCCA kernels */
  occa::kernel divergenceKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel gradientKernel;
  occa::kernel lowerPressureKernel;
  occa::kernel raisePressureKernel;
  occa::kernel stiffnessKernel;
  occa::kernel vecScaleKernel;
  occa::kernel vecScaledAddKernel;
  occa::kernel vecZeroKernel;
  occa::kernel weightedInnerProductKernel;
  occa::kernel globalWeightedInnerProductKernel;
  occa::kernel stokesOperatorKernel;
  occa::kernel multipleGlobalWeightedInnerProductsKernel;
  
  /* Scratch variables */
  dlong NblockV;          /* Used for reductions over the velocity DOFs. */
  dfloat *workV;
  occa::memory o_workV;

  dlong NblockP;          /* Used for reductions over the pressure DOFs. */
  dfloat *workP;
  occa::memory o_workP;

  // workspace for stokes solve
  int NsolveWorkspace;
  occa::memory *o_solveWorkspace; //

  int NpreconWorkspace;
  occa::memory *o_preconWorkspace; // 
  
} stokes_t;

stokes_t *stokesSetup(occa::properties &kernelInfoV, occa::properties &kernelInfoP, setupAide options);
void stokesSolveSetup(stokes_t *stokes, dfloat lambda, dfloat *eta, occa::properties &kernelInfoV, occa::properties &kernelInfoP);
void stokesSolve(stokes_t *stokes, dfloat lambda);
void stokesTimeDependentSolve(stokes_t *stokes, dfloat tfinal);
void stokesOperator(stokes_t *stokes, dfloat lambda,  stokesVec_t v, stokesVec_t Av);
void stokesPreconditioner(stokes_t *stokes, dfloat lambda, stokesVec_t v, stokesVec_t Mv);
void stokesPreconditionerSetup(stokes_t *stokes, dfloat lambda, occa::properties &kernelInfoV, occa::properties &kernelInfoP);

void stokesVecAllocate(stokes_t *stokes, stokesVec_t *v);
void stokesVecFree(stokes_t *stokes, stokesVec_t *v);
void stokesVecCopyHostToDevice(stokesVec_t v);
void stokesVecCopyDeviceToHost(stokesVec_t v);

void stokesVecCopy(stokes_t *stokes, stokesVec_t u, stokesVec_t v);
void stokesVecGatherScatter(stokes_t *stokes, stokesVec_t v);
void stokesVecUnmaskedGatherScatter(stokes_t *stokes, stokesVec_t v);
void stokesVecInnerProduct(stokes_t *stokes, stokesVec_t u, stokesVec_t v, dfloat *c);
void stokesVecScale(stokes_t *stokes, stokesVec_t v, dfloat c);
void stokesVecScaledAdd(stokes_t *stokes, dfloat a, stokesVec_t u, dfloat b, stokesVec_t v);
void stokesVecZero(stokes_t *stokes, stokesVec_t v);

void stokesOperatorPrint(stokes_t *stokes, dfloat lambda);
void stokesVecPrint(stokes_t *stokes, stokesVec_t v);

void stokesGetTestCase(stokes_t *stokes, stokesTestCase_t *testCase);

void stokesPlotVTU(stokes_t *stokes, char *fileName);

// occa memory versions

void stokesOperator(stokes_t *stokes, dfloat lambda, occa::memory &v, occa::memory &Av);
void stokesVecInnerProduct(stokes_t *stokes, occa::memory &u, occa::memory &v, dfloat *c);
void stokesVecCopy(stokes_t *stokes, occa::memory &u, occa::memory &v);

void stokesVecScaledAdd(stokes_t *stokes,
			dfloat a, occa::memory &u,
			dfloat b, occa::memory &v);

void stokesVecScale(stokes_t *stokes, occa::memory &v, dfloat c);

void stokesPreconditioner(stokes_t *stokes, dfloat lambda, occa::memory &v, occa::memory &Mv);

void stokesSolve(stokes_t *stokes, dfloat lambda, occa::memory &f, occa::memory &u);

#endif /* STOKES_H */
