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

int main(int argc, char **argv) {
  int              Nvelocity, Npressure, dim, elementType;
  string           fileName;
  mesh_t           *meshV, *meshP;
  stokes_t         *stokes;
  occa::properties kernelInfoV, kernelInfoP;

  // Start up MPI.
  MPI_Init(&argc, &argv);

  if (argc != 2) {
    printf("usage: ./stokesMain setupfile\n");
    MPI_Finalize();
    exit(-1);
  }

  // Load information from the setup file.
  setupAide options(argv[1]);
  options.getArgs("MESH FILE", fileName);
  options.getArgs("VELOCITY POLYNOMIAL DEGREE", Nvelocity);
  options.getArgs("PRESSURE POLYNOMIAL DEGREE", Npressure);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  // Setup meshes.
  //
  // TODO:  Is two separate meshes the right approach?
  meshV = meshSetup((char*)fileName.c_str(), Nvelocity, options);
  meshP = meshSetup((char*)fileName.c_str(), Npressure, options);

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
      meshOccaSetup2D(meshV, options, kernelInfoV);
      meshOccaSetup2D(meshP, options, kernelInfoP);
    } else {
      /* TODO:  Implement. */
      printf("ERROR:  Support for elements other than QUADRILATERALS not yet implemented.\n");
      MPI_Finalize();
      exit(-1);
    }
  } else if (dim == 3) {
    /* TODO:  Implement. */
    printf("ERROR:  Support for 3D elements not yet implemented.\n");
    MPI_Finalize();
    exit(-1);
  } else {
    printf("ERROR:  MESH DIMENSION must be 2 or 3.\n");
    MPI_Finalize();
    exit(-1);
  }

  stokes = stokesSetup(meshV, meshP, kernelInfoP, kernelInfoV, options);

  /* Solve. */
  
  /* Compute error (if applicable.) */

  /* Export solution. */
  
  /* Report runtime statistics. */

  // Shut down MPI.
  MPI_Finalize();

  return 0;
}
