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

static void stokesApplyStiffness(stokes_t *stokes, dfloat *v, dfloat *Av);
static void stokesApplyDivergenceX(stokes_t *stokes, dfloat *v, dfloat *Av);
static void stokesApplyDivergenceY(stokes_t *stokes, dfloat *v, dfloat *Av);
static void stokesApplyDivergenceXTranspose(stokes_t *stokes, dfloat *v, dfloat *Av);
static void stokesApplyDivergenceYTranspose(stokes_t *stokes, dfloat *v, dfloat *Av);

/* Applies the Stokes operator to the vector v, storing the result in Av.  Note
 * that v and Av cannot be the same.
 *
 * TODO:  This is woefully inefficient---giant wastes of both flops and
 * memory---but we're going to blow it all away once we get the thing working
 * anyway.
 *
 * TODO:  This only works for 2D quadrilaterals at the moment.
 */
void stokesOperator(stokes_t *stokes, stokesVec_t v, stokesVec_t Av)
{
  dfloat *Kvx, *Kvy, *Bxvp, *Byvp, *BxTvx, *ByTvy;

  Kvx = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  Kvy = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  Bxvp = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  Byvp = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  BxTvx = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));
  ByTvy = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));

  /* TODO:  Replace with calls to appropriate kernels (e.g., elliptic Ax kernel). */
  stokesApplyStiffness(stokes, v.x, Kvx);
  stokesApplyStiffness(stokes, v.y, Kvy);
  stokesApplyDivergenceX(stokes, v.p, Bxvp);
  stokesApplyDivergenceY(stokes, v.p, Byvp);
  stokesApplyDivergenceXTranspose(stokes, v.x, BxTvx);
  stokesApplyDivergenceYTranspose(stokes, v.y, ByTvy);

  /* Add up the contributions from the blocks. */
  for (int i = 0; i < stokes->NtotalV; i++) {
    Av.x[i] = 0;
    Av.y[i] = 0;
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    Av.p[i] = 0;
  }

  for (int i = 0; i < stokes->NtotalV; i++) {
    Av.x[i] += Kvx[i] - Bxvp[i];
    Av.y[i] += Kvy[i] - Byvp[i];
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    Av.p[i] = -BxTvx[i] - ByTvy[i];
  }

  free(Kvx);
  free(Kvy);
  free(Bxvp);
  free(Byvp);
  free(BxTvx);
  free(ByTvy);

  // Gather-scatter for C0 FEM.
  //
  // TODO:  Make a function for this.
  if (stokes->options.compareArgs("VELOCITY DISCRETIZATION", "CONTINUOUS")) {
    ogsGatherScatter(Av.x, ogsDfloat, ogsAdd, stokes->meshV->ogs);
    ogsGatherScatter(Av.y, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  }

  if (stokes->options.compareArgs("PRESSURE DISCRETIZATION", "CONTINUOUS")) {
    ogsGatherScatter(Av.p, ogsDfloat, ogsAdd, stokes->meshP->ogs);
  }

  return;
}

static void stokesApplyDivergenceXTranspose(stokes_t *stokes, dfloat *v, dfloat *Av)
{
  mesh_t *meshV, *meshP;
  dfloat *tmp0, *tmp1, *tmpITy0, *tmpITy1, *Av0, *Av1;

  meshV = stokes->meshV;
  meshP = stokes->meshP;

  tmp0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmpITy0 = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  tmpITy1 = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  Av0 = (dfloat*)calloc(meshP->Np, sizeof(dfloat));
  Av1 = (dfloat*)calloc(meshP->Np, sizeof(dfloat));

  for (int e = 0; e < meshP->Nelements; e++) {
    memset(tmp0, 0, meshV->Np*sizeof(dfloat));
    memset(tmp1, 0, meshV->Np*sizeof(dfloat));
    memset(tmpITy0, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(tmpITy1, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(Av0, 0, meshP->Np*sizeof(dfloat));
    memset(Av1, 0, meshP->Np*sizeof(dfloat));

    /* Multiply by I (X) D or D (X) I. */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          tmp0[i*meshV->Nq + j] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i*meshV->Nq + k];
          tmp1[i + j*meshV->Nq] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i + k*meshV->Nq];
        }
      }
    }

    /* Apply geometric factors and quadrature weights. */
    for (int i = 0; i < meshV->Np; i++) {
      tmp0[i] *= meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*RXID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
      tmp1[i] *= meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*SXID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
    }

    /* Apply transpose of interpolation operator. */

    /* y-action */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          tmpITy0[i*meshV->Nq + j] += meshP->interpRaise[i + k*meshP->Nq]*tmp0[j + k*meshV->Nq];
          tmpITy1[i*meshV->Nq + j] += meshP->interpRaise[i + k*meshP->Nq]*tmp1[j + k*meshV->Nq];
        }
      }
    }

    /* x-action */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshP->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          Av0[i*meshP->Nq + j] += meshP->interpRaise[j + k*meshP->Nq]*tmpITy0[i*meshV->Nq + k];
          Av1[i*meshP->Nq + j] += meshP->interpRaise[j + k*meshP->Nq]*tmpITy1[i*meshV->Nq + k];
        }
      }
    }

    /* Sum up to get the result. */
    for (int i = 0; i < meshP->Np; i++) {
      Av[e*meshP->Np + i] = Av0[i] + Av1[i];
    }
  }

  free(tmp0);
  free(tmp1);
  free(tmpITy0);
  free(tmpITy1);
  free(Av0);
  free(Av1);

  return;
}

static void stokesApplyDivergenceYTranspose(stokes_t *stokes, dfloat *v, dfloat *Av)
{
  mesh_t *meshV, *meshP;
  dfloat *tmp0, *tmp1, *tmpITy0, *tmpITy1, *Av0, *Av1;

  meshV = stokes->meshV;
  meshP = stokes->meshP;

  tmp0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmpITy0 = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  tmpITy1 = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  Av0 = (dfloat*)calloc(meshP->Np, sizeof(dfloat));
  Av1 = (dfloat*)calloc(meshP->Np, sizeof(dfloat));

  for (int e = 0; e < meshP->Nelements; e++) {
    memset(tmp0, 0, meshV->Np*sizeof(dfloat));
    memset(tmp1, 0, meshV->Np*sizeof(dfloat));
    memset(tmpITy0, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(tmpITy1, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(Av0, 0, meshP->Np*sizeof(dfloat));
    memset(Av1, 0, meshP->Np*sizeof(dfloat));

    /* Multiply by I (X) D or D (X) I. */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          tmp0[i*meshV->Nq + j] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i*meshV->Nq + k];
          tmp1[i + j*meshV->Nq] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i + k*meshV->Nq];
        }
      }
    }

    /* Apply geometric factors and quadrature weights. */
    for (int i = 0; i < meshV->Np; i++) {
      tmp0[i] *= meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*RYID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
      tmp1[i] *= meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*SYID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
    }

    /* Apply transpose of interpolation operator. */

    /* y-action */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          tmpITy0[i*meshV->Nq + j] += meshP->interpRaise[i + k*meshP->Nq]*tmp0[j + k*meshV->Nq];
          tmpITy1[i*meshV->Nq + j] += meshP->interpRaise[i + k*meshP->Nq]*tmp1[j + k*meshV->Nq];
        }
      }
    }

    /* x-action */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshP->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          Av0[i*meshP->Nq + j] += meshP->interpRaise[j + k*meshP->Nq]*tmpITy0[i*meshV->Nq + k];
          Av1[i*meshP->Nq + j] += meshP->interpRaise[j + k*meshP->Nq]*tmpITy1[i*meshV->Nq + k];
        }
      }
    }

    /* Sum up to get the result. */
    for (int i = 0; i < meshP->Np; i++) {
      Av[e*meshP->Np + i] = Av0[i] + Av1[i];
    }
  }

  free(tmp0);
  free(tmp1);
  free(tmpITy0);
  free(tmpITy1);
  free(Av0);
  free(Av1);

  return;
}

static void stokesApplyDivergenceX(stokes_t *stokes, dfloat *v, dfloat *Av)
{
  mesh_t *meshV, *meshP;
  dfloat *tmpIx, *tmpIxy, *tmp0, *tmp1, *Av0, *Av1;

  meshV = stokes->meshV;
  meshP = stokes->meshP;

  tmpIx = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  tmpIxy = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));

  for (int e = 0; e < meshP->Nelements; e++) {
    memset(tmpIx, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(tmpIxy, 0, meshV->Np*sizeof(dfloat));
    memset(tmp0, 0, meshV->Np*sizeof(dfloat));
    memset(tmp1, 0, meshV->Np*sizeof(dfloat));
    memset(Av0, 0, meshV->Np*sizeof(dfloat));
    memset(Av1, 0, meshV->Np*sizeof(dfloat));

    /* Interpolate to the higher-degree grid.
     *
     * TODO:  This assumes that the pressure space is one degree lower than the
     * velocity space.  In principle, we want to be more flexible than this.
     */

    /* x-interpolation */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshP->Nq; k++) {
          tmpIx[i*meshV->Nq + j] += meshP->interpRaise[j*meshP->Nq + k]*v[e*meshP->Np + i*meshP->Nq + k];
        }
      }
    }

    /* y-interpolation */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshP->Nq; k++) {
          tmpIxy[i*meshV->Nq + j] += meshP->interpRaise[i*meshP->Nq + k]*tmpIx[j + k*meshV->Nq];
        }
      }
    }

    /* Apply geometric factors and quadrature weights. */
    for (int i = 0; i < meshV->Np; i++) {
      tmp0[i] = tmpIxy[i]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*RXID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
      tmp1[i] = tmpIxy[i]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*SXID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
    }

    /* Multiply by I (X) D^T or D^T (X). */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          Av0[i*meshV->Nq + j] += meshV->D[k*meshV->Nq + j]*tmp0[i*meshV->Nq + k];
          Av1[i + j*meshV->Nq] += meshV->D[k*meshV->Nq + j]*tmp1[i + k*meshV->Nq];
        }
      }
    }

    /* Sum up to get the result. */
    for (int i = 0; i < meshV->Np; i++) {
      Av[e*meshV->Np + i] = Av0[i] + Av1[i];
    }

  }

  free(Av0);
  free(Av1);

  free(tmpIx);
  free(tmpIxy);
  free(tmp0);
  free(tmp1);

  return;
}

static void stokesApplyDivergenceY(stokes_t *stokes, dfloat *v, dfloat *Av)
{
  mesh_t *meshV, *meshP;
  dfloat *tmpIx, *tmpIxy, *tmp0, *tmp1, *Av0, *Av1;

  meshV = stokes->meshV;
  meshP = stokes->meshP;

  tmpIx = (dfloat*)calloc(meshV->Nq*meshP->Nq, sizeof(dfloat));
  tmpIxy = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av0 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av1 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));

  for (int e = 0; e < meshP->Nelements; e++) {
    memset(tmpIx, 0, meshV->Nq*meshP->Nq*sizeof(dfloat));
    memset(tmpIxy, 0, meshV->Np*sizeof(dfloat));
    memset(tmp0, 0, meshV->Np*sizeof(dfloat));
    memset(tmp1, 0, meshV->Np*sizeof(dfloat));
    memset(Av0, 0, meshV->Np*sizeof(dfloat));
    memset(Av1, 0, meshV->Np*sizeof(dfloat));

    /* Interpolate to the higher-degree grid.
     *
     * TODO:  This assumes that the pressure space is one degree lower than the
     * velocity space.  In principle, we want to be more flexible than this.
     */

    /* x-interpolation */
    for (int i = 0; i < meshP->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshP->Nq; k++) {
          tmpIx[i*meshV->Nq + j] += meshP->interpRaise[j*meshP->Nq + k]*v[e*meshP->Np + i*meshP->Nq + k];
        }
      }
    }

    /* y-interpolation */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshP->Nq; k++) {
          tmpIxy[i*meshV->Nq + j] += meshP->interpRaise[i*meshP->Nq + k]*tmpIx[j + k*meshV->Nq];
        }
      }
    }

    /* Apply geometric factors and quadrature weights. */
    for (int i = 0; i < meshV->Np; i++) {
      tmp0[i] = tmpIxy[i]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*RYID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
      tmp1[i] = tmpIxy[i]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*SYID]*meshV->vgeo[e*meshV->Nvgeo*meshV->Np + i + meshV->Np*JWID];
    }

    /* Multiply by I (X) D^T or D^T (X). */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          Av0[i*meshV->Nq + j] += meshV->D[k*meshV->Nq + j]*tmp0[i*meshV->Nq + k];
          Av1[i + j*meshV->Nq] += meshV->D[k*meshV->Nq + j]*tmp1[i + k*meshV->Nq];
        }
      }
    }

    /* Sum up to get the result. */
    for (int i = 0; i < meshV->Np; i++) {
      Av[e*meshV->Np + i] = Av0[i] + Av1[i];
    }

  }

  free(Av0);
  free(Av1);

  free(tmpIx);
  free(tmpIxy);
  free(tmp0);
  free(tmp1);

  return;
}

static void stokesApplyStiffness(stokes_t *stokes, dfloat *v, dfloat *Av)
{
  mesh_t *meshV, *meshP;
  dfloat *Av00, *Av01, *Av10, *Av11;
  dfloat *tmp00, *tmp01, *tmp10, *tmp11;

  meshV = stokes->meshV;
  meshP = stokes->meshP;

  Av00 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av01 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av10 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  Av11 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));

  tmp00 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp01 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp10 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));
  tmp11 = (dfloat*)calloc(meshV->Np, sizeof(dfloat));

  for (int e = 0; e < meshV->Nelements; e++) {
    memset(Av00, 0, meshV->Np*sizeof(dfloat));
    memset(Av01, 0, meshV->Np*sizeof(dfloat));
    memset(Av10, 0, meshV->Np*sizeof(dfloat));
    memset(Av11, 0, meshV->Np*sizeof(dfloat));

    memset(tmp00, 0, meshV->Np*sizeof(dfloat));
    memset(tmp01, 0, meshV->Np*sizeof(dfloat));
    memset(tmp10, 0, meshV->Np*sizeof(dfloat));
    memset(tmp11, 0, meshV->Np*sizeof(dfloat));

    /* Apply I (X) D or D (X) I as appropriate. */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          tmp00[i*meshV->Nq + j] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i*meshV->Nq + k];
          tmp01[i + j*meshV->Nq] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i + k*meshV->Nq];
          tmp10[i*meshV->Nq + j] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i*meshV->Nq + k];
          tmp11[i + j*meshV->Nq] += meshV->D[j*meshV->Nq + k]*v[e*meshV->Np + i + k*meshV->Nq];
        }
      }
    }

    /* Apply the geometric factors, the Jacobian, and the quadrature weights. */
    for (int i = 0; i < meshV->Np; i++) {
        tmp00[i] *= meshV->ggeo[e*meshV->Nggeo*meshV->Np + i + meshV->Np*G00ID];
        tmp01[i] *= meshV->ggeo[e*meshV->Nggeo*meshV->Np + i + meshV->Np*G01ID];
        tmp10[i] *= meshV->ggeo[e*meshV->Nggeo*meshV->Np + i + meshV->Np*G01ID];
        tmp11[i] *= meshV->ggeo[e*meshV->Nggeo*meshV->Np + i + meshV->Np*G11ID];
    }

    /* Apply I (X) D^T or D^T (X) I as appropriate. */
    for (int i = 0; i < meshV->Nq; i++) {
      for (int j = 0; j < meshV->Nq; j++) {
        for (int k = 0; k < meshV->Nq; k++) {
          Av00[i*meshV->Nq + j] += meshV->D[k*meshV->Nq + j]*tmp00[i*meshV->Nq + k];
          Av01[i*meshV->Nq + j] += meshV->D[k*meshV->Nq + j]*tmp01[i*meshV->Nq + k];
          Av11[i + j*meshV->Nq] += meshV->D[k*meshV->Nq + j]*tmp11[i + k*meshV->Nq];
          Av10[i + j*meshV->Nq] += meshV->D[k*meshV->Nq + j]*tmp10[i + k*meshV->Nq];
        }
      }
    }

    /* Sum up to get the result. */
    for (int i = 0; i < meshV->Np; i++) {
      Av[e*meshV->Np + i] = Av00[i] + Av01[i] + Av10[i] + Av11[i];
    }
  }

  free(Av00);
  free(Av01);
  free(Av10);
  free(Av11);

  free(tmp00);
  free(tmp01);
  free(tmp10);
  free(tmp11);

  return;
}
