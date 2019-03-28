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

static void stokesJacobiPreconditionerSetup(stokes_t *stokes);
static void stokesSchurComplementBlockDiagPreconditionerSetup(stokes_t *stokes, dfloat lambda, occa::properties &kernelInfoV, occa::properties & kernelInfoP);

static void stokesBuildLocalContinuousDiagQuad2D(stokes_t* stokes, dlong e, dfloat *diagA);

void stokesPreconditionerSetup(stokes_t *stokes, dfloat lambda, occa::properties &kernelInfoV, occa::properties &kernelInfoP)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokes->precon = NULL;
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    // TODO:  Update for lambda.
    if (lambda != 0.0) {
      printf("ERROR:  JACOBI preconditioner not implemented for nonzero LAMBDA.\n");
      exit(-1);
    }

    stokesJacobiPreconditionerSetup(stokes);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "SCHURCOMPLEMENTBLOCKDIAG")) {
    stokesSchurComplementBlockDiagPreconditionerSetup(stokes, lambda, kernelInfoV, kernelInfoP);
  } else {
    printf("ERROR:  Invalid value %s for [PRECONDITIONER] option.\n",
           stokes->options.getArgs("PRECONDITIONER").c_str());
    exit(-1);
  }

  return;
}

static void stokesJacobiPreconditionerSetup(stokes_t *stokes)
{
  int elementType;

  stokes->precon = new stokesPrecon_t();

  stokes->options.getArgs("JACOBI BOOST PARAMETER", stokes->precon->boost);
  if (stokes->precon->boost == 0.0) {
    printf("ERROR:  Jacobi boost must be nonzero.\n");
    exit(-1);
  }

  stokesVecAllocate(stokes, &stokes->precon->invDiagA);

  switch (stokes->elementType) {
  case QUADRILATERALS:
    if (stokes->meshV->dim == 2) {
      for (dlong e = 0; e < stokes->meshV->Nelements; e++) {
        /* TODO:  This is wasteful---should just build once and then copy. */
        stokesBuildLocalContinuousDiagQuad2D(stokes, e, stokes->precon->invDiagA.x + e*stokes->meshV->Np);
        stokesBuildLocalContinuousDiagQuad2D(stokes, e, stokes->precon->invDiagA.y + e*stokes->meshV->Np);
      }

      for (int i = 0; i < stokes->NtotalP; i++)
        stokes->precon->invDiagA.p[i] = stokes->precon->boost;
    } else {
      printf("ERROR:  Not implemented.\n");
      exit(-1);
    }
    break;
  default:
    printf("ERROR:  Not implemented.\n");
    exit(-1);
  }
  
  for (int i = 0; i < stokes->Ndof; i++)
    stokes->precon->invDiagA.v[i] = 1.0/stokes->precon->invDiagA.v[i];
  
  stokesVecCopyHostToDevice(stokes->precon->invDiagA);

  if (stokes->options.compareArgs("VELOCITY DISCRETIZATION", "CONTINUOUS")) {
    stokesVecGatherScatter(stokes, stokes->precon->invDiagA);
  }

  return;
}

static void stokesSchurComplementBlockDiagPreconditionerSetup(stokes_t *stokes, dfloat lambda, occa::properties &kernelInfoV, occa::properties &kernelInfoP)
{
  mesh_t     *meshV;
  elliptic_t *elliptic;
  setupAide  ellipticOptions;

  meshV = stokes->meshV;
  elliptic = new elliptic_t();
  stokes->precon = new stokesPrecon_t();

  /* Set up the elliptic sub-solver. */
  elliptic->mesh = meshV;
  elliptic->dim = elliptic->mesh->dim;
  elliptic->elementType = stokes->elementType;

  /* TODO:  Map these to the Stokes setup file. */

  ellipticOptions.setArgs("BASIS", "NODAL");
  ellipticOptions.setArgs("DISCRETIZATION", "CONTINUOUS");
  ellipticOptions.setArgs("DEBUG ENABLE OGS", "1");
  ellipticOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");
  ellipticOptions.setArgs("KRYLOV SOLVER", "PCG");

  if (stokes->options.compareArgs("VELOCITY BLOCK PRECONDITIONER", "MULTIGRID")) {
    ellipticOptions.setArgs("PRECONDITIONER", "MULTIGRID");
    ellipticOptions.setArgs("MULTIGRID COARSENING",           stokes->options.getArgs("VELOCITY MULTIGRID COARSENING"));
    ellipticOptions.setArgs("MULTIGRID SMOOTHER",             stokes->options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
    ellipticOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",     stokes->options.getArgs("VELOCITY MULTIGRID CHEBYSHEV DEGREE"));
    ellipticOptions.setArgs("PARALMOND AGGREGATION STRATEGY", stokes->options.getArgs("VELOCITY MULTIGRID PARALMOND AGGREGATION STRATEGY"));
    ellipticOptions.setArgs("PARALMOND CYCLE",                stokes->options.getArgs("VELOCITY MULTIGRID PARALMOND CYCLE"));
    ellipticOptions.setArgs("PARALMOND SMOOTHER",             stokes->options.getArgs("VELOCITY MULTIGRID PARALMOND SMOOTHER"));
    ellipticOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",     stokes->options.getArgs("VELOCITY MULTIGRID PARALMOND CHEBYSHEV DEGREE"));
  } else if (stokes->options.compareArgs("VELOCITY BLOCK PRECONDITIONER", "JACOBI")) {
    ellipticOptions.setArgs("PRECONDITIONER", "JACOBI");
  }

  ellipticOptions.setArgs("VERBOSE", "FALSE");
  ellipticOptions.setArgs("INTEGRATION TYPE", stokes->options.getArgs("INTEGRATION TYPE"));
  ellipticOptions.setArgs("ELLIPTIC INTEGRATION", stokes->options.getArgs("INTEGRATION TYPE"));
  elliptic->options = ellipticOptions;

  elliptic->BCType = stokes->BCType;

  elliptic->r = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  elliptic->o_r = elliptic->mesh->device.malloc(stokes->NtotalV*sizeof(dfloat), elliptic->r);
  elliptic->x = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  elliptic->o_x = elliptic->mesh->device.malloc(stokes->NtotalV*sizeof(dfloat), elliptic->x);

  /* TODO:  This allocates a whole lot of extra stuff---we may be able to share
   * some scratch arrays with the stokes_t.
   */
  ellipticSolveSetup(elliptic, lambda, kernelInfoV);

  if (stokes->options.compareArgs("VELOCITY BLOCK PRECONDITIONER", "MULTIGRID"))
    parAlmond::Report(elliptic->precon->parAlmond);

  stokes->precon->ellipticV = elliptic;

  if (stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MASSMATRIX") ||
      stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MULTIGRID")) {
    stokesVecAllocate(stokes, &stokes->precon->invMM);

    for (dlong e = 0; e < stokes->meshP->Nelements; e++) {
      for (int n = 0; n < stokes->meshP->Np; n++) {
        stokes->precon->invMM.p[e*stokes->meshP->Np + n] = stokes->meshP->ggeo[e*stokes->meshP->Np*stokes->meshP->Nggeo + GWJID*stokes->meshP->Np + n];
      }
    }

    // assemble mass matrix (suitable for C0 pressure)
    stokesVecCopyHostToDevice(stokes->precon->invMM);
    stokesVecGatherScatter(stokes, stokes->precon->invMM);
    stokesVecCopyDeviceToHost(stokes->precon->invMM);

    for (dlong e = 0; e < stokes->meshP->Nelements; e++) {
      for (int n = 0; n < stokes->meshP->Np; n++) {
        dfloat val = stokes->precon->invMM.p[e*stokes->meshP->Np + n];
        if (val) {
          stokes->precon->invMM.p[e*stokes->meshP->Np + n] = 1.0/val;
        } else if (val < 0) {
          printf("APA:  Got negative value on pressure mass matrix diagonal!\n");
          exit(-1);
        } else {
          printf("APA:  Got zero on pressure mass matrix diagonal!\n");
          exit(-1);
        }
      }
    }
    
    stokesVecCopyHostToDevice(stokes->precon->invMM);
  }

  if (stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MULTIGRID")) {
    setupAide ellipticOptionsP;

    stokes->precon->ellipticP = new elliptic_t();
    stokes->precon->ellipticP->mesh = stokes->meshP;
    stokes->precon->ellipticP->dim = stokes->precon->ellipticP->mesh->dim;
    stokes->precon->ellipticP->elementType = stokes->elementType;

    ellipticOptionsP.setArgs("BASIS", "NODAL");
    ellipticOptionsP.setArgs("DISCRETIZATION", "CONTINUOUS");
    ellipticOptionsP.setArgs("DEBUG ENABLE OGS", "1");
    ellipticOptionsP.setArgs("DEBUG ENABLE REDUCTIONS", "1");
    ellipticOptionsP.setArgs("KRYLOV SOLVER", "PCG");
    ellipticOptionsP.setArgs("PRECONDITIONER", "MULTIGRID");

    ellipticOptionsP.setArgs("MULTIGRID COARSENING",           stokes->options.getArgs("PRESSURE MULTIGRID COARSENING"));
    ellipticOptionsP.setArgs("MULTIGRID SMOOTHER",             stokes->options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
    ellipticOptionsP.setArgs("MULTIGRID CHEBYSHEV DEGREE",     stokes->options.getArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE"));
    ellipticOptionsP.setArgs("PARALMOND AGGREGATION STRATEGY", stokes->options.getArgs("PRESSURE MULTIGRID PARALMOND AGGREGATION STRATEGY"));
    ellipticOptionsP.setArgs("PARALMOND CYCLE",                stokes->options.getArgs("PRESSURE MULTIGRID PARALMOND CYCLE"));
    ellipticOptionsP.setArgs("PARALMOND SMOOTHER",             stokes->options.getArgs("PRESSURE MULTIGRID PARALMOND SMOOTHER"));
    ellipticOptionsP.setArgs("PARALMOND CHEBYSHEV DEGREE",     stokes->options.getArgs("PRESSURE MULTIGRID PARALMOND CHEBYSHEV DEGREE"));

    ellipticOptionsP.setArgs("VERBOSE", "FALSE");
    ellipticOptionsP.setArgs("INTEGRATION TYPE",     stokes->options.getArgs("INTEGRATION TYPE"));
    ellipticOptionsP.setArgs("ELLIPTIC INTEGRATION", stokes->options.getArgs("INTEGRATION TYPE"));

    stokes->precon->ellipticP->options = ellipticOptionsP;

    int BCType[3] = {0, 2, 2};
    stokes->precon->ellipticP->BCType = (int*)calloc(3, sizeof(int));
    memcpy(stokes->precon->ellipticP->BCType, BCType, 3*sizeof(int));

    stokes->precon->ellipticP->r = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));
    stokes->precon->ellipticP->o_r = stokes->precon->ellipticP->mesh->device.malloc(stokes->NtotalP*sizeof(dfloat), elliptic->r);
    stokes->precon->ellipticP->x = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));
    stokes->precon->ellipticP->o_x = elliptic->mesh->device.malloc(stokes->NtotalP*sizeof(dfloat), elliptic->x);

    /* TODO:  This allocates a whole lot of extra stuff---we may be able to share
     * some scratch arrays with the stokes_t.
     */
    ellipticSolveSetup(stokes->precon->ellipticP, 0.0, kernelInfoP);

    if (stokes->options.compareArgs("PRESSURE BLOCK PRECONDITIONER", "MULTIGRID"))
      parAlmond::Report(stokes->precon->ellipticP->precon->parAlmond);
  }

  return;
}

/*****************************************************************************/

/* TODO:  This was basically copied from the elliptic solver.
 *
 * TODO:  This may need modification for boundary conditions.
 */
static void stokesBuildLocalContinuousDiagQuad2D(stokes_t* stokes, dlong e, dfloat *diagA)
{
  mesh_t *mesh = stokes->meshV;

  for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      int iid = nx+ny*mesh->Nq;
      diagA[iid] = 0;

      for (int k=0;k<mesh->Nq;k++) {
        int id = k+ny*mesh->Nq;
        dfloat eta = stokes->eta[e*mesh->Np + id];
        dfloat Grr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];
        diagA[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq]*eta;
      }

      for (int k=0;k<mesh->Nq;k++) {
        int id = nx+k*mesh->Nq;
        dfloat eta = stokes->eta[e*mesh->Np + id];
        dfloat Gss = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];
        diagA[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq]*eta;
      }

      int id = nx+ny*mesh->Nq;
      dfloat eta = stokes->eta[e*mesh->Np + id];
      dfloat Grs = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
      diagA[iid] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq]*eta;
    }
  }

  return;
}
