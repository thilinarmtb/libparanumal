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

  stokes = new stokes_t();

  // Load information from the setup file.
  options.getArgs("MESH FILE", fileName);
  options.getArgs("VELOCITY POLYNOMIAL DEGREE", velocityN);
  options.getArgs("PRESSURE POLYNOMIAL DEGREE", pressureN);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  stokes->options = options;

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
      printf("ERROR:  Support for elements other than QUADRILATERALS not yet implemented.\n");
      MPI_Finalize();
      exit(-1);
    }
  } else if (dim == 3) {
    printf("ERROR:  Support for 3D elements not yet implemented.\n");
    MPI_Finalize();
    exit(-1);
  } else {
    printf("ERROR:  MESH DIMENSION must be 2 or 3.\n");
    MPI_Finalize();
    exit(-1);
  }

  stokesSolveSetup(stokes, kernelInfoV, kernelInfoP);
  stokesSetupRHS(stokes);

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
      dfloat x, y, J;

      ind = e*stokes->meshV->Np + i;
      x = stokes->meshV->x[ind];
      y = stokes->meshV->y[ind];

      stokes->f.x[ind] = cos(M_PI*x)*sin(M_PI*y) - 24.0*y*pow(1.0 - x*x, 3.0)*(5.0*y*y - 3.0) - 36.0*y*(1.0 - x*x)*(5.0*x*x - 1.0)*pow(1.0 - y*y, 2.0);
      stokes->f.y[ind] = sin(M_PI*x)*cos(M_PI*y) + 24.0*x*pow(1.0 - y*y, 3.0)*(5.0*x*x - 3.0) + 36.0*x*(1.0 - y*y)*(5.0*y*y - 1.0)*pow(1.0 - x*x, 2.0);

      // NB:  We have to incorporate the Jacobian factor because meshApplyElementMatrix() assumes it.
      //
      // TODO:  This may assume the use of quadrilateral elements.
      J = stokes->meshV->vgeo[stokes->meshV->Np*(e*stokes->meshV->Nvgeo + JID) + i];
      stokes->f.x[ind] *= J;
      stokes->f.y[ind] *= J;
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
