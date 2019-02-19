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

/* Applies the Stokes operator to the vector v, storing the result in Av.  Note
 * that v and Av cannot be the same.
 *
 * TODO:  This only works for Quad2D at the moment.
 */
void stokesOperator(stokes_t *stokes, stokesVec_t v, stokesVec_t Av)
{
  /* TODO:  We re-allocate these scratch variables every time we call the
   * operator, but they're going to go away eventually, so we don't care.  keep
   * re-allocating them.
   */
  //occa::memory o_interpRaise = stokes->meshV->device.malloc(stokes->meshP->Nq*stokes->meshV->Nq*sizeof(dfloat), stokes->meshP->interpRaise);
  //occa::memory o_pRaised = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat));

  occa::memory o_Pp = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat));

  stokesVecZero(stokes, Av);

  /* NB:  These kernels MUST be called in this order, as they modify Av incrementally.
   *
   * TODO:  Fuse these into one big Stokes Ax kernel?
   */

  stokes->stiffnessKernel(stokes->meshV->Nelements,
                          stokes->meshV->o_ggeo,
                          stokes->meshV->o_Dmatrices,
                          stokes->o_eta,
                          v.o_x,
                          Av.o_x);

  stokes->stiffnessKernel(stokes->meshV->Nelements,
                          stokes->meshV->o_ggeo,
                          stokes->meshV->o_Dmatrices,
                          stokes->o_eta,
                          v.o_y,
                          Av.o_y);

  if (stokes->meshV->dim == 3) {
    stokes->stiffnessKernel(stokes->meshV->Nelements,
                            stokes->meshV->o_ggeo,
                            stokes->meshV->o_Dmatrices,
                            stokes->o_eta,
                            v.o_z,
                            Av.o_z);
  }

  /*
  stokes->raisePressureKernel(stokes->meshV->Nelements,
                              o_interpRaise,
                              v.o_p,
                              o_pRaised);
  */

  stokes->pressureProjectKernel(stokes->meshV->Nelements,
                                stokes->o_P,
                                v.o_p,
                                o_Pp);

  stokes->gradientKernel(stokes->meshV->Nelements,
                         stokes->NtotalV,
                         stokes->meshV->o_Dmatrices,
                         stokes->meshV->o_vgeo,
                         o_Pp,
                         Av.o_v);

  stokes->divergenceKernel(stokes->meshV->Nelements,
                           stokes->NtotalV,
                           stokes->meshV->o_Dmatrices,
                           stokes->meshV->o_vgeo,
                           v.o_v,
                           o_Pp);

  stokes->pressureProjectTransKernel(stokes->meshV->Nelements,
                                     stokes->o_P,
                                     o_Pp,
                                     Av.o_p);

  /*
  stokes->lowerPressureKernel(stokes->meshV->Nelements,
                              o_interpRaise,
                              o_pRaised,
                              Av.o_p);
  */

#if 0
  /* Rank-boost for all-Neumann problem.
   *
   * TODO:  Need to make this conditional on an "allNeumann" flag like in the
   * elliptic solver.
   *
   * TODO:  Do this on-device.
   *
   * TODO:  There may be a better way to do this using projections---see the
   * experimental branch.
   *
   */

  stokesVecCopyDeviceToHost(v);
  stokesVecCopyDeviceToHost(Av);

  dfloat sx = 0.0, sy = 0.0, sz = 0.0;
  for (int i = 0; i < stokes->NtotalV; i++) {
    sx += v.x[i];
    sy += v.y[i];
    if (stokes->meshV->dim == 3)
      sz += v.z[i];
  }

  sx /= stokes->NtotalV;
  sy /= stokes->NtotalV;
  if (stokes->meshV->dim == 3)
    sz /= stokes->NtotalV;

  for (int i = 0; i < stokes->NtotalV; i++) {
    Av.x[i] += sx;
    Av.y[i] += sy;
    if (stokes->meshV->dim == 3)
      Av.z[i] += sz;
  }

  stokesVecCopyHostToDevice(Av);
#endif

  /* Gather-scatter for C0 FEM. */
  if (stokes->options.compareArgs("VELOCITY DISCRETIZATION", "CONTINUOUS"))
    stokesVecGatherScatter(stokes, Av);

  // TODO:  Make a function for this.
  //
  // TODO:  We only need to do this for C0 FEM.
  if (stokes->Nmasked) {
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_x);
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_y);
    if (stokes->meshV->dim == 3)
      stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_z);
  }

  o_Pp.free();
  //o_pRaised.free();
  //o_interpRaise.free();
  return;
}

/*****************************************************************************/

void stokesOperatorPrint(stokes_t *stokes)
{
  stokesVec_t v, Av;

  stokesVecAllocate(stokes, &v);
  stokesVecAllocate(stokes, &Av);

  for (int i = 0; i < stokes->Ndof; i++) {
    v.v[i] = 1.0;
    stokesVecCopyHostToDevice(v);
    stokesOperator(stokes, v, Av);
    stokesVecCopyDeviceToHost(Av);
    stokesVecPrint(stokes, Av);
    v.v[i] = 0.0;
  }

  stokesVecFree(stokes, &v);
  stokesVecFree(stokes, &Av);
}
