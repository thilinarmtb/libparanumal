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
  stokes_t         *stokes;
  occa::properties kernelInfoV, kernelInfoP;

  // Start up MPI.
  MPI_Init(&argc, &argv);

  if (argc != 2) {
    printf("usage: ./stokesMain setupfile\n");
    MPI_Finalize();
    exit(-1);
  }

  setupAide options(argv[1]);

  stokes = stokesSetup(kernelInfoP, kernelInfoV, options);
  stokesSolve(stokes);

  /* Compute error (if applicable.) */

  /* Export solution. */
  
  /* Report runtime statistics. */

  // Shut down MPI.
  MPI_Finalize();

  return 0;
}
