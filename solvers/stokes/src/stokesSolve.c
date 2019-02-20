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

void stokesSolve(stokes_t *stokes)
{
  if (stokes->options.compareArgs("KRYLOV SOLVER", "MINRES")) {
    stokesSolveMINRES(stokes);
  } else {
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

    // HERE ==+>
#if 1
    dfloat alpha =  -(del/gam);
    dfloat beta  = (i==0) ? 0: -(gam/gamp);
    stokesUpdateMINRES(stokes, -a2, -a3, alpha, beta, z, w_old, w, r_old, r, p);
#else
    stokesVecScaledAdd(stokes, -a2, w, 1.0, z);                /* z = z - a2*w - a3*w_old  */
    stokesVecScaledAdd(stokes, -a3, w_old, 1.0, z);
    stokesVecCopy(stokes, w, w_old);                           /* w_old = w                */
    stokesVecCopy(stokes, z, w);                               /* w = z                    */
    stokesVecCopy(stokes, r, z);                               /* z = r                    */
    stokesVecScaledAdd(stokes, 1.0, p, -(del/gam), r);         /* r = p - (del/gam)*r      */
    if (i > 0)
      stokesVecScaledAdd(stokes, -(gam/gamp), r_old, 1.0, r);  /* r = r - (gam/gamp)*r_old */
    stokesVecCopy(stokes, z, r_old);                           /* r_old = z                */
    // TO HERE <+=====
#endif

    stokesPreconditioner(stokes, r, z);                        /* z = M\r                  */
    gamp = gam;
    stokesVecInnerProduct(stokes, z, r, &gam);                 /* gam = sqrt(r'*z)         */
    gam = sqrt(gam);
    a1 = sqrt(a0*a0 + gam*gam);
    cp = c;
    c  = a0/a1;
    sp = s;
    s  = gam/a1;

    // FROM HERE
    stokesVecScale(stokes, w, 1.0/a1);                         /* w = w/a1                 */
    stokesVecScaledAdd(stokes, c*eta, w, 1.0, u);              /* u = u + c*eta*w          */
    /// TO HERE
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
