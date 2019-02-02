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

static void stokesAllocateVec(stokes_t *stokes, stokesVec_t *v);
static void stokesFreeVec(stokes_t *stokes, stokesVec_t *v);

static void stokesVecScale(stokes_t *stokes, stokesVec_t v, dfloat c);
static void stokesVecScaledAdd(stokes_t *stokes, dfloat a, stokesVec_t y, dfloat b, stokesVec_t v);
static void stokesVecInnerProduct(stokes_t *stokes, stokesVec_t u, stokesVec_t v, dfloat *c);
static void stokesVecCopy(stokes_t *stokes, stokesVec_t u, stokesVec_t v);

static void stokesVecPrint(stokes_t *stokes, stokesVec_t v)
{
  for (int i = 0; i < stokes->NtotalV; i++)
    printf("% .15e\n", v.x[i]);
  for (int i = 0; i < stokes->NtotalV; i++)
    printf("% .15e\n", v.y[i]);
  for (int i = 0; i < stokes->NtotalP; i++)
    printf("% .15e\n", v.p[i]);
}

void stokesSolve(stokes_t *stokes)
{
  stokesVec_t u, p, z, r, r_old, w, w_old;
  dfloat      a0, a1, a2, a3, del, gam, gamp, c, cp, s, sp, eta;

  u = stokes->u;

  /* Allocate work vectors. */
  stokesAllocateVec(stokes, &p);
  stokesAllocateVec(stokes, &z);
  stokesAllocateVec(stokes, &r);
  stokesAllocateVec(stokes, &r_old);
  stokesAllocateVec(stokes, &w);
  stokesAllocateVec(stokes, &w_old);

  stokesOperator(stokes, u, r);                                /* r = f - Au               */
  stokesVecScaledAdd(stokes, 1.0, stokes->f, -1.0, r);
  stokesPreconditioner(stokes, r, z);                          /* z = M\r                  */
  stokesVecInnerProduct(stokes, z, r, &gam);                   /* gam = sqrt(r'*z)         */
  gamp = 0.0;
  gam  = sqrt(gam);                           
  eta  = gam;
  sp   = 0.0;
  s    = 0.0;
  cp   = 1.0;
  c    = 1.0;

  /* TODO:  Get a proper convergence criterion. */
  for (int i = 0; i < 100; i++) {
    stokesVecScale(stokes, z, 1.0/gam);                        /* z = z/gam                */
    stokesOperator(stokes, z, p);                              /* p = Az                   */
    stokesVecInnerProduct(stokes, p, z, &del);                 /* del = z'*p               */
    a0 = c*del - cp*s*gam;
    a2 = s*del + cp*c*gam;
    a3 = sp*gam;
    stokesVecScaledAdd(stokes, -a2, w, 1.0, z);                /* z = z - a2*w - a3*w_old  */
    stokesVecScaledAdd(stokes, -a3, w_old, 1.0, z);
    stokesVecCopy(stokes, w, w_old);                           /* w_old = w                */
    stokesVecCopy(stokes, z, w);                               /* w = z                    */
    stokesVecCopy(stokes, r, z);                               /* z = r                    */
    stokesVecScaledAdd(stokes, 1.0, p, -(del/gam), r);         /* r = p - (del/gam)*r      */
    if (i > 0)
      stokesVecScaledAdd(stokes, -(gam/gamp), r_old, 1.0, r);  /* r = r - (gam/gamp)*r_old */
    stokesVecCopy(stokes, z, r_old);                           /* r_old = z                */
    stokesPreconditioner(stokes, r, z);                        /* z = M\r                  */
    gamp = gam;
    stokesVecInnerProduct(stokes, z, r, &gam);                 /* gam = sqrt(r'*z)         */
    gam = sqrt(gam);
    a1 = sqrt(a0*a0 + gam*gam);
    cp = c;
    c  = a0/a1;
    sp = s;
    s  = gam/a1;
    stokesVecScale(stokes, w, 1.0/a1);                         /* w = w/a1                 */
    stokesVecScaledAdd(stokes, c*eta, w, 1.0, u);              /* u = u + c*eta*w          */
    eta = -s*eta;

    printf("gam = % .15e\n", gam);
    if (fabs(gam) < 1.0e-10)
      break;
  }

  /* Free work vectors. */
  stokesFreeVec(stokes, &p);
  stokesFreeVec(stokes, &z);
  stokesFreeVec(stokes, &r);
  stokesFreeVec(stokes, &r_old);
  stokesFreeVec(stokes, &w);
  stokesFreeVec(stokes, &w_old);

  return;
}

static void stokesAllocateVec(stokes_t *stokes, stokesVec_t *v)
{
  v->x = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  v->y = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  if (stokes->meshV->dim == 3)
    v->z = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  else
    v->z = NULL;
  v->p = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));

  return;
}

static void stokesFreeVec(stokes_t *stokes, stokesVec_t *v)
{
  free(v->x);
  free(v->y);
  if (stokes->meshV->dim == 0)
    free(v->z);
  free(v->p);

  return;
}

// TODO:  Only works for 2D.
static void stokesVecScale(stokes_t *stokes, stokesVec_t v, dfloat c)
{
  for (int i = 0; i < stokes->NtotalV; i++) {
    v.x[i] *= c;
    v.y[i] *= c;
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    v.p[i] *= c;
  }

  return;
}

// Computes v <-- au + bv.
//
// TODO:  Only works for 2D.
static void stokesVecScaledAdd(stokes_t *stokes, dfloat a, stokesVec_t u, dfloat b, stokesVec_t v)
{
  for (int i = 0; i < stokes->NtotalV; i++) {
    v.x[i] = a*u.x[i] + b*v.x[i];
    v.y[i] = a*u.y[i] + b*v.y[i];
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    v.p[i] = a*u.p[i] + b*v.p[i];
  }
}

// Computes c = v'*u
//
// TODO:  Only works for 2D.
static void stokesVecInnerProduct(stokes_t *stokes, stokesVec_t u, stokesVec_t v, dfloat *c)
{
  *c = 0.0;
  for (int i = 0; i < stokes->NtotalV; i++) {
    *c += u.x[i]*v.x[i] + u.y[i]*v.y[i];
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    *c += u.p[i]*v.p[i];
  }
}

// Copies v <-- u.
//
// TODO:  Only works for 2D.
static void stokesVecCopy(stokes_t *stokes, stokesVec_t u, stokesVec_t v)
{
  for (int i = 0; i < stokes->NtotalV; i++) {
    v.x[i] = u.x[i];
    v.y[i] = u.y[i];
  }

  for (int i = 0; i < stokes->NtotalP; i++) {
    v.p[i] = u.p[i];
  }
}
