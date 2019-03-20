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

static void stokesSolveMINRES(stokes_t *stokes);
static void stokesSolveDQGMRES(stokes_t *stokes);

void stokesSolve(stokes_t *stokes)
{
  if (stokes->options.compareArgs("KRYLOV SOLVER", "MINRES")) {
    stokesSolveMINRES(stokes);
  } else if (stokes->options.compareArgs("KRYLOV SOLVER", "DQGMRES")) {
    stokesSolveDQGMRES(stokes);
  }else {
    printf("ERROR:  Invalid value %s for [KRYLOV SOLVER] option.",
           stokes->options.getArgs("KRYLOV SOLVER").c_str());
    exit(-1);
  }
}

static void stokesSolveMINRES(stokes_t *stokes)
{
  stokesVec_t u, p, z, r, r_old, w, w_old;
  dfloat      a0, a1, a2, a3, del, gam, gamp, c, cp, s, sp, eta;
  dfloat      maxiter, tol;
  int         verbose;

  stokes->options.getArgs("KRYLOV SOLVER ITERATION LIMIT", maxiter);
  stokes->options.getArgs("KRYLOV SOLVER TOLERANCE",       tol);

  verbose = 0;
  if (stokes->options.compareArgs("VERBOSE", "TRUE"))
    verbose = 1;

  u = stokes->u;

  /* Allocate work vectors.
   *
   * TODO:  These vectors reside entirely on the device---no need to waste
   * memory for their host counterparts.
   *
   * TODO:  Also:  should these be owned by the stokes_t so we don't have to
   * re-allocate every time we want to solve?
   */

  stokesVecAllocate(stokes, &p);
  stokesVecAllocate(stokes, &z);
  stokesVecAllocate(stokes, &r);
  stokesVecAllocate(stokes, &r_old);
  stokesVecAllocate(stokes, &w);
  stokesVecAllocate(stokes, &w_old);

  stokesOperator(stokes, u, r);                                /* r = f - Au               */
  stokesVecScaledAdd(stokes, 1.0, stokes->f, -1.0, r);
  stokesPreconditioner(stokes, r, z);                          /* z = M\r                  */
  stokesVecInnerProduct(stokes, z, r, &gam);                   /* gam = sqrt(r'*z)         */
  gamp = 0.0;
  if (gam < 0) {
    printf("BAD:  gam < 0.\n");
    exit(-1);
  }
  gam  = sqrt(gam);
  eta  = gam;
  sp   = 0.0;
  s    = 0.0;
  cp   = 1.0;
  c    = 1.0;

  /* Adjust the tolerance to account for small initial residual norms. */
  tol = mymax(tol*fabs(eta), tol);
  if (verbose)
    printf("MINRES:  initial eta = % .15e, target %.15e\n", eta, tol);

  /* MINRES iteration loop. */
  for (int i = 0; i < maxiter; i++) {
    if (verbose)
      printf("MINRES:  it % 3d  eta = % .15e\n", i, eta);
    if (fabs(eta) < tol) {
      if (verbose)
        printf("MINRES converged in %d iterations (eta = % .15e).\n", i, eta);
      break;
    }

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
  }

  /* Free work vectors. */
  stokesVecFree(stokes, &p);
  stokesVecFree(stokes, &z);
  stokesVecFree(stokes, &r);
  stokesVecFree(stokes, &r_old);
  stokesVecFree(stokes, &w);
  stokesVecFree(stokes, &w_old);

  return;
}

static void stokesSolveDQGMRES(stokes_t *stokes)
{
  stokesVec_t u, f, e, w, p1, p2, p3, v1, v2, v3, tmp;
  dfloat      g1, g2, c1, s1, c2, s2, c3, s3, h0, h1, h2, h3, a;
  dfloat      maxiter, tol;
  int         verbose;

  stokes->options.getArgs("KRYLOV SOLVER ITERATION LIMIT", maxiter);
  stokes->options.getArgs("KRYLOV SOLVER TOLERANCE",       tol);

  verbose = 0;
  if (stokes->options.compareArgs("VERBOSE", "TRUE"))
    verbose = 1;

  u = stokes->u;
  f = stokes->f;

  /* Allocate work vectors.
   *
   * TODO:  These vectors reside entirely on the device---no need to waste
   * memory for their host counterparts.
   *
   * TODO:  Put these in the stokes_t so we don't re-allocate them every time
   * we want to solve.
   *
   * TODO:  Can we eliminate two of these so that we use the same storage as MINRES?
   */
  stokesVecAllocate(stokes, &e);
  stokesVecAllocate(stokes, &w);
  stokesVecAllocate(stokes, &v1);
  stokesVecAllocate(stokes, &v2);
  stokesVecAllocate(stokes, &v3);
  stokesVecAllocate(stokes, &p1);
  stokesVecAllocate(stokes, &p2);
  stokesVecAllocate(stokes, &p3);

  stokesOperator(stokes, u, v2);                          /* v2 = f - A*u0    (initial residual) */
  stokesVecScaledAdd(stokes, 1.0, f, -1.0, v2);
  stokesVecInnerProduct(stokes, v2, v2, &g1);             /* g1 = norm(v2)                       */
  g1 = sqrt(g1);
  stokesVecScale(stokes, v2, 1.0/g1);                     /* v2 = v2/g1                          */

  /* Adjust the tolerance to account for small initial residual norm. */
  tol = mymax(tol*fabs(g1), tol);
  if (verbose)
    printf("DQGMRES:  initial gamma = % .15e, target %.15e\n", g1, tol);

  /* DQGMRES iteration loop. */
  for (int i = 0; i < maxiter; i++) {
    if (verbose)
      printf("DQGMRES:  it % 3d  gamma = % .15e\n", i, g1);
    if (fabs(g1) < tol) {
      if (verbose)
        printf("DQGMRES converged in %d iterations (gamma = % .15e).\n", i, g1);
      break;
    }

    /* Compute initial choice for the next vector. */
    stokesPreconditioner(stokes, v2, w);                  /* w = M\v2                            */
    stokesOperator(stokes, w, v3);                        /* v3 = Aw                             */

    /* Orthogonalize (incompletely). */
    if (i > 0) {
      stokesVecInnerProduct(stokes, v3, v1, &h1);         /* h1 = v1'*v3                         */
      stokesVecScaledAdd(stokes, -h1, v1, 1.0, v3);       /* v3 = v3 - h1*v1                     */
    }
    stokesVecInnerProduct(stokes, v3, v2, &h2);           /* h2 = v3'v2                          */
    stokesVecScaledAdd(stokes, -h2, v2, 1.0, v3);         /* v3 = v3 - h2*v2                     */

    /* Normalize. */
    stokesVecInnerProduct(stokes, v3, v3, &h3);           /* h3 = norm(v3)                       */
    h3 = sqrt(h3);
    stokesVecScale(stokes, v3, 1.0/h3);                   /* v3 = v3/h3                          */

    /* Apply previous rotations to the new column of the Hessenberg matrix. */
    if (i > 1) {
      h0 = s1*h1;
      h1 = c1*h1;
    }

    if (i > 0) {
      dfloat h1n, h2n;
      h1n =  c2*h1 + s2*h2;
      h2n = -s2*h1 + c2*h2;
      h1 = h1n;
      h2 = h2n;
    }

    /* Compute rotation to eliminate the new subdiagonal entry. */
    a = hypot(h2, h3);
    if (h3 != 0.0) {
      c3 = h2/a;
      s3 = h3/a;
    } else {
      c3 = 1.0;
      s3 = 0.0;
    }

    /* Apply new rotation to the Hessenberg matrix. */
    h2 = a;
    h3 = 0.0;

    /* Apply new rotation to the right-hand side. */
    g2 = -s3*g1;
    g1 =  c3*g1;

    /* Update the solution.
     *
     * TODO:  This can be made more efficient.
     */
    stokesVecCopy(stokes, w, p3);                         /* p3 = w                              */
    if (i > 1)
      stokesVecScaledAdd(stokes, -h0, p1, 1.0, p3);       /* p3 = p3 - h0*p1                     */
    if (i > 0)
      stokesVecScaledAdd(stokes, -h1, p2, 1.0, p3);       /* p3 = p3 - h1*p2                     */
    stokesVecScale(stokes, p3, 1.0/h2);                   /* p3 = p3/h2                          */

    stokesVecScaledAdd(stokes, g1, p3, 1.0, e);           /* e = e + g1*p3                       */

    /* Rotate the variables.
     *
     * NB:  Vector assignments copy *pointers*, not values.
     */
    tmp = v1;
    v1 = v2;
    v2 = v3;
    v3 = tmp;

    tmp = p1;
    p1 = p2;
    p2 = p3;
    p3 = tmp;

    c1 = c2;
    s1 = s2;

    c2 = c3;
    s2 = s3;

    g1 = g2;
  }

  /* Undo the right preconditioning. */
  stokesVecScaledAdd(stokes, 1.0, e, 1.0, u);

  /* Free work vectors. */
  stokesVecFree(stokes, &e);
  stokesVecFree(stokes, &w);
  stokesVecFree(stokes, &v1);
  stokesVecFree(stokes, &v2);
  stokesVecFree(stokes, &v3);
  stokesVecFree(stokes, &p1);
  stokesVecFree(stokes, &p2);
  stokesVecFree(stokes, &p3);

  return;
}
