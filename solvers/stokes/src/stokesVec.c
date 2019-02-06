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

void stokesVecAllocate(stokes_t *stokes, stokesVec_t *v)
{
  v->v = (dfloat*)calloc(stokes->Ndof, sizeof(dfloat));

  v->x = v->v;
  v->y = v->v + stokes->NtotalV;
  if (stokes->meshV->dim == 2) {
    v->z = NULL;
    v->p = v->v + 2*stokes->NtotalV;
  } else if (stokes->meshV->dim == 3) {
    v->z = v->v + 2*stokes->NtotalV;
    v->p = v->v + 3*stokes->NtotalV;
  }

  return;
}

/* Copies v <-- u.
 *
 * TODO:  Replace with OCCA kernel.
 */
void stokesVecCopy(stokes_t *stokes, stokesVec_t u, stokesVec_t v)
{
  for (int i = 0; i < stokes->Ndof; i++)
    v.v[i] = u.v[i];

  return;
}

void stokesVecGatherScatter(stokes_t *stokes, stokesVec_t v)
{
  ogsGatherScatter(v.x, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  ogsGatherScatter(v.y, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  if (stokes->meshV->dim == 3)
    ogsGatherScatter(v.z, ogsDfloat, ogsAdd, stokes->meshV->ogs);
  ogsGatherScatter(v.p, ogsDfloat, ogsAdd, stokes->meshP->ogs);

  return;
}

void stokesVecFree(stokes_t *stokes, stokesVec_t *v)
{
  free(v->v);

  v->v = NULL;
  v->x = NULL;
  v->y = NULL;
  v->z = NULL;
  v->p = NULL;

  return;
}

/* Computes v <-- c*v for vector v, scalar c.
 *
 * TODO:  Replace with OCCA kernel.
 */
void stokesVecScale(stokes_t *stokes, stokesVec_t v, dfloat c)
{
  for (int i = 0; i < stokes->Ndof; i++)
    v.v[i] *= c;

  return;
}

/* Computes v <-- au + bv for vectors u, v and scalars a, b.
 *
 * TODO:  Replace with OCCA kernel.
 */
void stokesVecScaledAdd(stokes_t *stokes, dfloat a, stokesVec_t u, dfloat b, stokesVec_t v)
{
  for (int i = 0; i < stokes->Ndof; i++)
    v.v[i] = a*u.v[i] + b*v.v[i];

  return;
}

/* Computes c = v'*u
 *
 * TODO:  The inverse degree weighting is only applicable for C0 FEM.
 *
 * TODO:  Replace with OCCA kernel.
 */
void stokesVecInnerProduct(stokes_t *stokes, stokesVec_t u, stokesVec_t v, dfloat *c)
{
  *c = 0.0;

  /* Easier to go component-by-component because of the weights. */
  if (stokes->meshV->dim == 2) {
    for (int i = 0; i < stokes->NtotalV; i++)
      *c += (u.x[i]*v.x[i] + u.y[i]*v.y[i])*stokes->meshV->ogs->invDegree[i];
  } else if (stokes->meshV->dim == 3) {
    for (int i = 0; i < stokes->NtotalV; i++)
      *c += (u.x[i]*v.x[i] + u.y[i]*v.y[i] + u.z[i]*v.z[i])*stokes->meshV->ogs->invDegree[i];
  }

  for (int i = 0; i < stokes->NtotalP; i++)
    *c += u.p[i]*v.p[i]*stokes->meshP->ogs->invDegree[i];

  return;
}

/* Sets v <-- 0. */
void stokesVecZero(stokes_t *stokes, stokesVec_t v)
{
  memset(v.v, 0, stokes->Ndof*sizeof(dfloat));
  return;
}

/*****************************************************************************/

void stokesVecPrint(stokes_t *stokes, stokesVec_t v)
{
  for (int i = 0; i < stokes->Ndof; i++)
    printf("% .20e\n", v.v[i]);

  return;
}
