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
  const dfloat zero = 0.0;
  const dfloat one  = 1.0;
  const dfloat mone = -1.0;
  const dfloat tau  = 0.0;

  /* TODO:  We re-allocate these scratch variables every time we call the
   * operator, but they're going to go away eventually, so we don't care.  keep
   * re-allocating them.
   */
  occa::memory o_pProjected = stokes->mesh->device.malloc(stokes->Ntotal*sizeof(dfloat));

  stokesVecZero(stokes, Av);

  /* NB:  These kernels MUST be called in this order, as they modify Av incrementally.
   *
   * TODO:  Fuse these into one big Stokes Ax kernel?
   */

  if(!stokes->options.compareArgs("STIFFNESS INTEGRATION TYPE", "CUBATURE")){
    printf("APA:  Old world (solve/stiffness).\n");
    stokes->stiffnessKernel(stokes->mesh->Nelements,
			    stokes->mesh->o_ggeo,
			    stokes->mesh->o_Dmatrices,
			    stokes->o_eta,
			    v.o_x,
			    Av.o_x);
    
    stokes->stiffnessKernel(stokes->mesh->Nelements,
			    stokes->mesh->o_ggeo,
			    stokes->mesh->o_Dmatrices,
			    stokes->o_eta,
			    v.o_y,
			    Av.o_y);
    
    if (stokes->mesh->dim == 3) {
      stokes->stiffnessKernel(stokes->mesh->Nelements,
			      stokes->mesh->o_ggeo,
			      stokes->mesh->o_Dmatrices,
			      stokes->o_eta,
			      v.o_z,
			      Av.o_z);
    }
  }
  else{
    
    printf("APA:  New world (solve/stiffness).\n");
    stokes->stiffnessKernel(stokes->mesh->Nelements,
			    stokes->mesh->o_cubggeo,
			    stokes->o_cubD,
			    stokes->o_cubInterp,
			    stokes->o_cubEta,
			    v.o_x,
			    Av.o_x);
    
    stokes->stiffnessKernel(stokes->mesh->Nelements,
			    stokes->mesh->o_cubggeo,
			    stokes->o_cubD,
			    stokes->o_cubInterp,
			    stokes->o_cubEta,
                          v.o_y,
			    Av.o_y);
    
    if (stokes->mesh->dim == 3) {
      stokes->stiffnessKernel(stokes->mesh->Nelements,
			      stokes->mesh->o_cubggeo,
			      stokes->o_cubD,
			      stokes->o_cubInterp,
			      stokes->o_cubEta,
                            v.o_z,
			      Av.o_z);
    }
  }


  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
                                  one,
                                  mone,
                                  stokes->o_uP,
                                  stokes->o_vP,
                                  v.o_p,
                                  o_pProjected);

  if(!stokes->options.compareArgs("STIFFNESS INTEGRATION TYPE", "CUBATURE")){
    printf("APA:  Old world (solve/grad-div).\n");
    stokes->gradientKernel(stokes->mesh->Nelements,
			   stokes->Ntotal,
			   stokes->mesh->o_Dmatrices,
			   stokes->mesh->o_vgeo,
			   o_pProjected,
			   Av.o_v);
    
    stokes->divergenceKernel(stokes->mesh->Nelements,
			     stokes->Ntotal,
			     stokes->mesh->o_Dmatrices,
			     stokes->mesh->o_vgeo,
			     v.o_v,
			     o_pProjected);
  }else{
    
    printf("APA:  New world (solve/grad-div).\n");
    stokes->gradientKernel(stokes->mesh->Nelements,
			   stokes->Ntotal,
			   stokes->mesh->o_cubvgeo,
			   stokes->o_cubD,
			   stokes->o_cubInterp,
			   stokes->o_cubInterp,
			   o_pProjected,
			   Av.o_v);
    
    stokes->divergenceKernel(stokes->mesh->Nelements,
			     stokes->Ntotal,
			     stokes->mesh->o_cubvgeo,
			     stokes->o_cubD,
			     stokes->o_cubInterp,
			     stokes->o_cubInterp,
			     v.o_v,
			     o_pProjected);
  }
  
  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
                                  one,
                                  mone,
                                  stokes->o_vP,
                                  stokes->o_uP,
                                  o_pProjected,
                                  Av.o_p);

#if 0
  // Pressure stabilization block.
  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
                                  zero,
                                  one,
                                  stokes->o_uP,
                                  stokes->o_vP,
                                  v.o_p,
                                  o_pProjected);

  stokes->rankOneProjectionKernel(stokes->mesh->Nelements,
                                  zero,
                                  tau,
                                  stokes->o_vP,
                                  stokes->o_uP,
                                  o_pProjected,
                                  o_pProjected);

  stokes->vecScaledAddKernel(stokes->Ntotal,
                             one,
                             o_pProjected,
                             one,
                             Av.o_p);
#endif

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
  for (int i = 0; i < stokes->Ntotal; i++) {
    sx += v.x[i];
    sy += v.y[i];
    if (stokes->mesh->dim == 3)
      sz += v.z[i];
  }

  sx /= stokes->Ntotal;
  sy /= stokes->Ntotal;
  if (stokes->mesh->dim == 3)
    sz /= stokes->Ntotal;

  for (int i = 0; i < stokes->Ntotal; i++) {
    Av.x[i] += sx;
    Av.y[i] += sy;
    if (stokes->mesh->dim == 3)
      Av.z[i] += sz;
  }

  stokesVecCopyHostToDevice(Av);
#endif

  /* Gather-scatter for C0 FEM. */
  if (stokes->options.compareArgs("VELOCITY DISCRETIZATION", "CONTINUOUS")){
    stokesVecGatherScatter(stokes, Av);
#if 0
    ogsGatherScatter(Av.o_x, ogsDfloat, ogsAdd, stokes->mesh->ogs);
    ogsGatherScatter(Av.o_y, ogsDfloat, ogsAdd, stokes->mesh->ogs);
    if (stokes->mesh->dim == 3)
      ogsGatherScatter(Av.o_z, ogsDfloat, ogsAdd, stokes->mesh->ogs);
#endif
  }

  // TODO:  Make a function for this.
  //
  // TODO:  We only need to do this for C0 FEM.
  if (stokes->Nmasked) {
    stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_x);
    stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_y);
    if (stokes->mesh->dim == 3)
      stokes->mesh->maskKernel(stokes->Nmasked, stokes->o_maskIds, Av.o_z);
  }

  o_pProjected.free();
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
