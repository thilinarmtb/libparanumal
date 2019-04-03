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

  mesh_t *meshV = stokes->meshV;
  
  stokesVec_t tmp, tmp2;

  occa::memory o_pRaised = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat));

  stokesVecAllocate(stokes, &tmp);
  stokesVecAllocate(stokes, &tmp2);

  stokes->userBoundaryConditionsKernel(meshV->Nelements, stokes->NtotalV, t, stokes->o_mapB, meshV->o_x, meshV->o_y, meshV->o_z, tmp.o_v);
  
  if (stokes->options.compareArgs("INTEGRATION TYPE", "GLL")) {
    stokes->raisePressureKernel(stokes->meshV->Nelements,
                                stokes->o_interpRaise,
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
                                stokes->o_interpRaise,
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

  // Gather-scatter for C0 FEM.
  stokesVecUnmaskedGatherScatter(stokes, stokes->f);

  // TODO:  We only need to do this for C0 FEM.
  if (stokes->Nmasked) {
    stokes->velocityMaskKernel(stokes->Nmasked, stokes->NtotalV, stokes->o_maskIds, stokes->f.o_v);
    stokes->velocityMaskKernel(stokes->Nmasked, stokes->NtotalV, stokes->o_maskIds, stokes->u.o_v);
  }

  return;
}
