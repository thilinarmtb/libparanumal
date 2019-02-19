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

static void stokesSetupRHS(stokes_t *stokes);
static void stokesRHSAddBC(stokes_t *stokes);

static void stokesTestForcingFunctionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy);
static void stokesTestForcingFunctionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy);
static void stokesTestForcingFunctionDirichletQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy);
static void stokesTestForcingFunctionLeakyCavityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy);
static void stokesTestForcingFunctionConstantViscosityHex3D(dfloat x, dfloat y, dfloat z, dfloat *fx, dfloat *fy, dfloat *fz);

stokes_t *stokesSetup(occa::properties &kernelInfoV, occa::properties &kernelInfoP, setupAide options)
{
  int      velocityN, pressureN, dim, elementType;
  int      velocityNtotal, pressureNtotal;
  string   fileName, bcHeaderFileName;
  stokes_t *stokes;
  dfloat   *eta;

  stokes = new stokes_t();

  // Load information from the setup file.
  options.getArgs("MESH FILE", fileName);
  options.getArgs("DATA FILE", bcHeaderFileName);
  options.getArgs("VELOCITY POLYNOMIAL DEGREE", velocityN);
  options.getArgs("PRESSURE POLYNOMIAL DEGREE", pressureN);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  stokes->options = options;
  stokes->elementType = elementType;

  // Setup meshes.
  //
  // TODO:  Is two separate meshes the right approach?
  //
  // TODO:  This results in duplicate device objects and instruction streams.
  stokes->meshV = meshSetup((char*)fileName.c_str(), velocityN, options);
  stokes->meshP = meshSetup((char*)fileName.c_str(), pressureN, options);

  // OCCA setup.
  kernelInfoV["defines"].asObject();
  kernelInfoV["includes"].asArray();
  kernelInfoV["header"].asArray();
  kernelInfoV["flags"].asObject();

  kernelInfoP["defines"].asObject();
  kernelInfoP["includes"].asArray();
  kernelInfoP["header"].asArray();
  kernelInfoP["flags"].asObject();

  kernelInfoV["includes"] += bcHeaderFileName.c_str();

  if (dim == 2) {
    if (elementType == QUADRILATERALS) {
      meshOccaSetup2D(stokes->meshV, options, kernelInfoV);
      meshOccaSetup2D(stokes->meshP, options, kernelInfoP);
    } else {
      printf("ERROR:  Support for 2D elements other than QUADRILATERALS not yet implemented.\n");
      MPI_Finalize();
      exit(-1);
    }
  } else if (dim == 3) {
    if (elementType == HEXAHEDRA) {
      meshOccaSetup3D(stokes->meshV, options, kernelInfoV);
      meshOccaSetup3D(stokes->meshP, options, kernelInfoP);
    } else {
      printf("ERROR:  Support for 3D elements other than HEXAHEDRA not yet implemented.\n");
      MPI_Finalize();
      exit(-1);
    }
  } else {
    printf("ERROR:  MESH DIMENSION must be 2 or 3.\n");
    MPI_Finalize();
    exit(-1);
  }

  /* TODO:  Think about where this should be set. */
  eta = (dfloat*)calloc(stokes->meshV->Nelements*stokes->meshV->Np, sizeof(dfloat));
  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->Np; i++) {
      int    ind;
      dfloat x, y, z;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];
      z = stokes->meshV->y[ind];

      if (dim == 2) {
        eta[ind] = 1.0;
        //eta[ind] = 2.0 + sinh(x*y);
      } else if (dim == 3) {
        eta[ind] = 1.0;
      }
    }
  }

  // Set up the physical-to-mathematial BC map.
  int BCType[3] = {0, 1, 2};
  stokes->BCType = (int*)calloc(3, sizeof(int));
  memcpy(stokes->BCType, BCType, 3*sizeof(int));

  stokesSolveSetup(stokes, eta, kernelInfoV, kernelInfoP);
  stokesSetupRHS(stokes);

  free(eta);
  return stokes;
}

static void stokesSetupRHS(stokes_t *stokes)
{
  int dim;

  stokes->options.getArgs("MESH DIMENSION", dim);

  // Initialize right-hand side with the forcing term.
  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->Np; i++) {
      int    ind;
      dfloat x, y, z, J;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];
      z = stokes->meshV->z[ind];

      if (dim == 2) {
        //stokesTestForcingFunctionConstantViscosityQuad2D(x, y, stokes->f.x + ind, stokes->f.y + ind);
        //stokesTestForcingFunctionVariableViscosityQuad2D(x, y, stokes->f.x + ind, stokes->f.y + ind);
        stokesTestForcingFunctionDirichletQuad2D(x, y, stokes->f.x + ind, stokes->f.y + ind);
      } else if (dim == 3) {
        stokesTestForcingFunctionConstantViscosityHex3D(x, y, z, stokes->f.x + ind, stokes->f.y + ind, stokes->f.z + ind);
      }

      // NB:  We have to incorporate the Jacobian factor because meshApplyElementMatrix() assumes it.
      //
      // TODO:  This may assume the use of quadrilateral/hexahedral elements.
      J = stokes->meshV->vgeo[stokes->meshV->Np*(e*stokes->meshV->Nvgeo + JID) + i];
      stokes->f.x[ind] *= J;
      stokes->f.y[ind] *= J;
      if (dim == 3)
        stokes->f.z[ind] *= J;
    }
  }

  // Multiply by mass matrix to get the true right-hand side.
  meshApplyElementMatrix(stokes->meshV, stokes->meshV->MM, stokes->f.x, stokes->f.x);
  meshApplyElementMatrix(stokes->meshV, stokes->meshV->MM, stokes->f.y, stokes->f.y);
  if (dim == 3)
    meshApplyElementMatrix(stokes->meshV, stokes->meshV->MM, stokes->f.z, stokes->f.z);

  // Move RHS to the device.
  stokesVecCopyHostToDevice(stokes->f);

  // Apply the boundary conditions.
  stokesRHSAddBC(stokes);

  // Gather-scatter for C0 FEM.
  stokesVecUnmaskedGatherScatter(stokes, stokes->f);

  printf("f vector before mask:\n");
  stokesVecCopyDeviceToHost(stokes->f);
  stokesVecPrint(stokes, stokes->f);

  // TODO:  Make a function for this.
  //
  // TODO:  We only need to do this for C0 FEM.
  if (stokes->Nmasked) {
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_x);
    stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_y);
    if (stokes->meshV->dim == 3)
      stokes->meshV->maskKernel(stokes->Nmasked, stokes->o_maskIds, stokes->f.o_z);
  }

  printf("f vector after mask:\n");
  stokesVecCopyDeviceToHost(stokes->f);
  stokesVecPrint(stokes, stokes->f);

  return;
}

// TODO:  Write a kernel for this.
static void stokesRHSAddBC(stokes_t *stokes)
{
  stokesVec_t tmp;

  occa::memory o_interpRaise = stokes->meshV->device.malloc(stokes->meshP->Nq*stokes->meshV->Nq*sizeof(dfloat), stokes->meshP->interpRaise);
  occa::memory o_pRaised = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat));

  stokesVecAllocate(stokes, &tmp);

  for (int e = 0; e < stokes->meshV->Nelements; e++) {
    for (int i = 0; i < stokes->meshV->Np; i++) {
      int    ind;
      dfloat x, y, z;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];

      if (stokes->mapB[ind] == 1) {
        tmp.x[ind] = cos(y);
        tmp.y[ind] = sin(x);
      }

      /*
      if (stokes->mapB[ind] == 1) {
        if (fabs(y - 1.0) < 1.0e-12) {
          tmp.x[ind] = 1.0;
          tmp.y[ind] = 0.0;
        } else {
          tmp.x[ind] = 0.0;
          tmp.y[ind] = 0.0;
        }
      }
      */

    }
  }

  stokesVecCopyHostToDevice(tmp);
  printf("tmp:\n");
  stokesVecPrint(stokes, tmp);

  //stokesVecCopy(stokes, stokes->f, tmp);

  stokesVecZero(stokes, stokes->u);

  stokes->stiffnessKernel(stokes->meshV->Nelements,
                          stokes->meshV->o_ggeo,
                          stokes->meshV->o_Dmatrices,
                          stokes->o_eta,
                          tmp.o_x,
                          stokes->u.o_x);

  stokes->stiffnessKernel(stokes->meshV->Nelements,
                          stokes->meshV->o_ggeo,
                          stokes->meshV->o_Dmatrices,
                          stokes->o_eta,
                          tmp.o_y,
                          stokes->u.o_y);

  if (stokes->meshV->dim == 3) {
    stokes->stiffnessKernel(stokes->meshV->Nelements,
                            stokes->meshV->o_ggeo,
                            stokes->meshV->o_Dmatrices,
                            stokes->o_eta,
                            tmp.o_z,
                            stokes->u.o_z);
  }

  stokes->raisePressureKernel(stokes->meshV->Nelements,
                              o_interpRaise,
                              tmp.o_p,
                              o_pRaised);

  stokes->gradientKernel(stokes->meshV->Nelements,
                         stokes->NtotalV,
                         stokes->meshV->o_Dmatrices,
                         stokes->meshV->o_vgeo,
                         o_pRaised,
                         stokes->u.o_v);

  stokes->divergenceKernel(stokes->meshV->Nelements,
                           stokes->NtotalV,
                           stokes->meshV->o_Dmatrices,
                           stokes->meshV->o_vgeo,
                           tmp.o_v,
                           o_pRaised);

  stokes->lowerPressureKernel(stokes->meshV->Nelements,
                              o_interpRaise,
                              o_pRaised,
                              stokes->u.o_p);


  stokesVecCopyDeviceToHost(stokes->u);
  stokesVecCopyDeviceToHost(stokes->f);

  printf("u before scaled add:\n");
  stokesVecPrint(stokes, stokes->u);
  printf("f before scaled add:\n");
  stokesVecPrint(stokes, stokes->f);


  stokesVecScaledAdd(stokes, -1.0, stokes->u, 1.0, stokes->f);


  stokesVecCopyDeviceToHost(stokes->f);
  printf("f after scaled add:\n");
  stokesVecPrint(stokes, stokes->f);

  stokesVecZero(stokes, stokes->u);

  stokesVecFree(stokes, &tmp);
  o_pRaised.free();
  o_interpRaise.free();
  return;
}


/*****************************************************************************/

static void stokesTestForcingFunctionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = cos(M_PI*x)*sin(M_PI*y)/M_PI - 24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) - 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0);
  *fy = sin(M_PI*x)*cos(M_PI*y)/M_PI + 24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0);
  return;
}

static void stokesTestForcingFunctionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = cos(M_PI*x)*sin(M_PI*y)/M_PI - (2.0 + sinh(x*y))*(24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) + 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0)) + 36.0*x*pow(1.0 - x*x, 2.0)*pow(1.0 - y*y, 2.0)*y*y*cosh(x*y) - 6.0*pow(1.0 - x*x, 3.0)*(1.0 - 5.0*y*y)*(1.0 - y*y)*x*cosh(x*y);
  *fy = sin(M_PI*x)*cos(M_PI*y)/M_PI + (2.0 + sinh(x*y))*(24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0)) + 6.0*pow(1.0 - y*y, 3.0)*(1.0 - 5.0*x*x)*(1.0 - x*x)*y*cosh(x*y) - 36.0*y*pow(1.0 - y*y, 2.0)*pow(1.0 - x*x, 2.0)*x*x*cosh(x*y);
  return;
}

static void stokesTestForcingFunctionDirichletQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = 1.0 + cos(y);
  *fy = 1.0 + sin(x);
  return;
}

static void stokesTestForcingFunctionLeakyCavityQuad2D(dfloat x, dfloat y, dfloat *fx, dfloat *fy)
{
  *fx = 0.0;
  *fy = 0.0;
  return;
}

static void stokesTestForcingFunctionConstantViscosityHex3D(dfloat x, dfloat y, dfloat z, dfloat *fx, dfloat *fy, dfloat *fz)
{
  *fx = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 48.0*z*z*z - 72.0*z*(1.0 - z*z);
  *fy = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z) + 48.0*x*x*x - 72.0*x*(1.0 - x*x);
  *fz = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z) + 48.0*y*y*y - 72.0*y*(1.0 - y*y);
  return;
}
