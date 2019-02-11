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

static void stokesBuildLocalContinuousDiagQuad2D(stokes_t* stokes, dlong e, dfloat *diagA);

void stokesPreconditionerSetup(stokes_t *stokes)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokes->precon = NULL;
    return;
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    stokesJacobiPreconditionerSetup(stokes);
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
