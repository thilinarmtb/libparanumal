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

#include <stokes.h>

static void stokesTimeDependentSolveBackwardEuler(stokes_t *stokes, dfloat tfinal);
static void stokesRHSAddBC(stokes_t *stokes, dfloat t, dfloat lambda);

void stokesTimeDependentSolve(stokes_t *stokes, dfloat tfinal)
{
  if (stokes->options.compareArgs("TIME INTEGRATOR", "BACKWARDEULER")) {
    stokesTimeDependentSolveBackwardEuler(stokes, tfinal);
  } else {
    printf("ERROR:  Invalid value %s for TIME INTEGRATOR option.\n",
           stokes->options.getArgs("TIME INTEGRATOR").c_str());
    MPI_Finalize();
    exit(-1);
  }

  return;
}

static void stokesTimeDependentSolveBackwardEuler(stokes_t *stokes, dfloat tfinal)
{
  dfloat t,   dt;
  stokesVec_t tmp;

  mesh_t *meshV = stokes->meshV;
  
  stokesVecAllocate(stokes, &tmp);

  stokes->options.getArgs("TIME STEP", dt);

  t = 0.0;
  while (t < tfinal) {
    /* Update the time. */
    t += dt;

    /* Set current RHS from current solution scaled by dt plus the forcing
     * function.
     */

    stokes->userForcingKernel(meshV->Nelements, stokes->NtotalV, t, 
			      meshV->o_x, meshV->o_y, meshV->o_z,
			      stokes->f.o_v);
    

    dfloat dtinv = 1./dt;
    stokes->prepareRhsKernel(meshV->Nelements, stokes->NtotalV, t, dtinv,
			     meshV->o_x, meshV->o_y, meshV->o_z,
			     meshV->o_vgeo, stokes->u.o_v, stokes->f.o_v);
    
#if 0
    stokesVecZero(stokes, stokes->f);

    stokesVecCopyDeviceToHost(stokes->f);
    stokesVecCopyDeviceToHost(stokes->u);

    for (int e = 0; e < meshV->Nelements; e++) {
      for (int i = 0; i < meshV->Np; i++) {
        int    ind;
        dfloat x, y, z, JW;

        ind = e*meshV->Np + i;
        x = meshV->x[ind];
        y = meshV->y[ind];
        z = meshV->z[ind];

        if (meshV->dim == 2)
          stokes->testCase->tdForcingFn2D(x, y, t, stokes->f.x + ind, stokes->f.y + ind);
        else if (meshV->dim == 3)
          stokes->testCase->tdForcingFn3D(x, y, z, t, stokes->f.x + ind, stokes->f.y + ind, stokes->f.z + ind);

        stokes->f.x[ind] += stokes->u.x[ind]/dt;
        stokes->f.y[ind] += stokes->u.y[ind]/dt;
        if (meshV->dim == 3)
          stokes->f.z[ind] += stokes->u.z[ind]/dt;

        /* Add the Jacobian factors (because meshApplyElementMatrix() demands it). */
        JW = meshV->vgeo[meshV->Np*(e*meshV->Nvgeo + JWID) + i];
        stokes->f.x[ind] *= JW;
        stokes->f.y[ind] *= JW;
        if (meshV->dim == 3)
          stokes->f.z[ind] *= JW;
      }
    }

    stokesVecCopyHostToDevice(stokes->f);
#endif
    
    /* Apply the boundary conditions. */
    stokesRHSAddBC(stokes, t, 1.0/dt);
    
    /* Solve Stokes system for new solution. */
    stokesSolve(stokes, 1.0/dt, stokes->f.o_v, stokes->u.o_v);

    printf("t = % .15e done\n", t);
  }

  stokesVecFree(stokes, &tmp);
  return;
}

// TODO:  This was copied from stokesSetup.c; need to put this into a kernel.
static void stokesRHSAddBC(stokes_t *stokes, dfloat t, dfloat lambda)
{
  stokesVec_t tmp, tmp2;

  occa::memory o_interpRaise = stokes->meshV->device.malloc(stokes->meshP->Nq*stokes->meshV->Nq*sizeof(dfloat), stokes->meshP->interpRaise);
  occa::memory o_pRaised = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat));

  stokesVecAllocate(stokes, &tmp);
  stokesVecAllocate(stokes, &tmp2);

  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->Np; i++) {
      int    ind;
      dfloat x, y, z, p;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];
      z = stokes->meshV->z[ind];

      /* TODO:  Handle boundary data when we don't know the exact solution. */
      if (stokes->mapB[ind] == 1) {
        if (stokes->meshV->dim == 2)
          stokes->testCase->tdSolFn2D(x, y, t, tmp.x + ind, tmp.y + ind, &p);
        else if (stokes->meshV->dim == 3)
          stokes->testCase->tdSolFn3D(x, y, z, t, tmp.x + ind, tmp.y + ind, tmp.z + ind, &p);
      }
    }
  }

  stokesVecCopyHostToDevice(tmp);

  if (stokes->options.compareArgs("INTEGRATION TYPE", "GLL")) {
    stokes->raisePressureKernel(stokes->meshV->Nelements,
                                o_interpRaise,
                                tmp.o_p,
                                o_pRaised);

    stokes->stokesOperatorKernel(stokes->meshV->Nelements,
                                 stokes->NtotalV,
                                 stokes->meshV->o_vgeo,
                                 stokes->meshV->o_Dmatrices,
                                 lambda,
                                 stokes->o_eta,
                                 tmp.o_v,
                                 o_pRaised,
                                 tmp2.o_v);

    stokes->lowerPressureKernel(stokes->meshV->Nelements,
                                o_interpRaise,
                                o_pRaised,
                                tmp2.o_p);
  } else if (stokes->options.compareArgs("INTEGRATION TYPE", "CUBATURE")) {
    stokes->stiffnessKernel(stokes->meshV->Nelements,
                            stokes->meshV->o_cubggeo,
                            stokes->o_cubD,
                            stokes->o_cubInterpV,
                            stokes->o_cubEta,
                            tmp.o_x,
                            tmp2.o_x);

    stokes->stiffnessKernel(stokes->meshV->Nelements,
                            stokes->meshV->o_cubggeo,
                            stokes->o_cubD,
                            stokes->o_cubInterpV,
                            stokes->o_cubEta,
                            tmp.o_y,
                            tmp2.o_y);

    if (stokes->meshV->dim == 3) {
    stokes->stiffnessKernel(stokes->meshV->Nelements,
                            stokes->meshV->o_cubggeo,
                            stokes->o_cubD,
                            stokes->o_cubInterpV,
                            stokes->o_cubEta,
                            tmp.o_z,
                            tmp2.o_z);
    }

    stokes->gradientKernel(stokes->meshV->Nelements,
			                     stokes->NtotalV,
			                     stokes->meshV->o_cubvgeo,
			                     stokes->o_cubD,
			                     stokes->o_cubInterpV,
			                     stokes->o_cubInterpP,
			                     tmp.o_p,
			                     tmp2.o_v);

    stokes->divergenceKernel(stokes->meshV->Nelements,
			                       stokes->NtotalV,
			                       stokes->meshV->o_cubvgeo,
			                       stokes->o_cubD,
			                       stokes->o_cubInterpV,
			                       stokes->o_cubInterpP,
			                       tmp.o_v,
			                       tmp2.o_p);
  }

  stokesVecScaledAdd(stokes, -1.0, tmp2, 1.0, stokes->f);

  stokesVecFree(stokes, &tmp);
  stokesVecFree(stokes, &tmp2);
  o_pRaised.free();
  o_interpRaise.free();

  // Gather-scatter for C0 FEM.
  stokesVecUnmaskedGatherScatter(stokes, stokes->f);

  // TODO:  Make a function for this.
  //
  // TODO:  We only need to do this for C0 FEM.
  if (stokes->Nmasked) {
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_x);
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_y);
    if (stokes->meshV->dim == 3)
      stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_z);
  }

  if (stokes->Nmasked) {
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->u.o_x);
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->u.o_y);
    if (stokes->meshV->dim == 3)
      stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->u.o_z);
  }

  return;
}
