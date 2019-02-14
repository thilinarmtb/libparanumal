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

stokes_t *stokesSetup(occa::properties &kernelInfoV, occa::properties &kernelInfoP, setupAide options)
{
  int      velocityN, pressureN, dim, elementType;
  int      velocityNtotal, pressureNtotal;
  string   fileName;
  stokes_t *stokes;
  dfloat   *eta;

  stokes = new stokes_t();

  // Load information from the setup file.
  options.getArgs("MESH FILE", fileName);
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

      //eta[ind] = 2.0 + sinh(x*y);
      eta[ind] = 1.0;
    }
  }

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

      // 2D square constant viscosity test case.
      //stokes->eta[ind] = 1.0;
      //stokes->f.x[ind] = cos(M_PI*x)*sin(M_PI*y)/M_PI - 24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) - 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0);
      //stokes->f.y[ind] = sin(M_PI*x)*cos(M_PI*y)/M_PI + 24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0);

      // 2D square variable viscosity test case.
      //stokes->eta[ind] = 2.0 + sinh(x*y);
      //stokes->f.x[ind] = cos(M_PI*x)*sin(M_PI*y)/M_PI - (2.0 + sinh(x*y))*(24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) + 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0)) + 36.0*x*pow(1.0 - x*x, 2.0)*pow(1.0 - y*y, 2.0)*y*y*cosh(x*y) - 6.0*pow(1.0 - x*x, 3.0)*(1.0 - 5.0*y*y)*(1.0 - y*y)*x*cosh(x*y);
      //stokes->f.y[ind] = sin(M_PI*x)*cos(M_PI*y)/M_PI + (2.0 + sinh(x*y))*(24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0)) + 6.0*pow(1.0 - y*y, 3.0)*(1.0 - 5.0*x*x)*(1.0 - x*x)*y*cosh(x*y) - 36.0*y*pow(1.0 - y*y, 2.0)*pow(1.0 - x*x, 2.0)*x*x*cosh(x*y);

      // 3D cube constant viscosity test case.
      //stokes->eta[ind] = 1.0;
      stokes->f.x[ind] = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 48.0*z*z*z - 72.0*z*(1.0 - z*z);
      stokes->f.y[ind] = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z) + 48.0*x*x*x - 72.0*x*(1.0 - x*x);
      stokes->f.z[ind] = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z) + 48.0*y*y*y - 72.0*y*(1.0 - y*y);

      // NB:  We have to incorporate the Jacobian factor because meshApplyElementMatrix() assumes it.
      //
      // TODO:  This may assume the use of quadrilateral/hexahedral elements.
      J = stokes->meshV->vgeo[stokes->meshV->Np*(e*stokes->meshV->Nvgeo + JID) + i];
      stokes->f.x[ind] *= J;
      stokes->f.y[ind] *= J;
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

  // Gather-scatter for C0 FEM.
  stokesVecGatherScatter(stokes, stokes->f);

  return;
}
