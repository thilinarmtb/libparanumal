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

#include <time.h>
#include "shell.h"

int main(int argc, char **argv)
{
  // Unbuffer stdout (can help with debugging).
  setbuf(stdout, NULL);

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./insMain setupfile\n");
    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  options.getArgs("MESH FILE", fileName);
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);
  
  // set up mesh
  mesh_t *mesh;
  if(dim==2){
    switch(elementType){   
    case TRIANGLES:
      mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
    case QUADRILATERALS:     
      mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
    }
  }
  else{
    dfloat radius = 1;

    switch(elementType){   
    case TRIANGLES:
      options.getArgs("SPHERE RADIUS", radius);
      mesh = meshSetupTri3D((char*)fileName.c_str(), N, radius); break;
    case QUADRILATERALS:
      options.getArgs("SPHERE RADIUS", radius);
      mesh = meshSetupQuad3D((char*)fileName.c_str(), N, radius); break;
    case TETRAHEDRA:
      mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
    case HEXAHEDRA:
      mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
    }
  }

  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  shell_t *shell = shellSetup(mesh, lambda, kernelInfo, options);

  // TODO:  If we just use mesh->device here instead of branching on the [SHELL
  // SOLVER] option, IPDG with MULTIGRID for SEM spits out NaNs.  Why?????
  occa::streamTag solveStart;
  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    solveStart = shell->mesh->device.tagStream();
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    solveStart = shell->elliptic->mesh->device.tagStream();
  }

  shellSolve(shell, options);

  occa::streamTag solveEnd;
  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    solveEnd = shell->mesh->device.tagStream();
    shell->times.solve = shell->mesh->device.timeBetween(solveStart, solveEnd);
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    solveEnd = shell->elliptic->mesh->device.tagStream();
    shell->times.solve = shell->elliptic->mesh->device.timeBetween(solveStart, solveEnd);
  }

  shell->times.total = shell->times.setup.total + shell->times.solve;

  if (!options.compareArgs("OUTPUT FILE NAME", "NONE")) {
    shellPlotVTU3D(shell, options.getArgs("OUTPUT FILE NAME").c_str(), 0);
  }

  dfloat errH1, errL2;
  shellErrorHex3D(shell, shell->q3D, &errH1, &errL2);

  printf("------------------------\n");
  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    printf("N\t\t%d\n",             shell->mesh->N);
    printf("Np\t\t%d\n",            shell->mesh->Np);
    printf("Nelements\t%d\n",       shell->mesh->Nelements);
    printf("Ntotal\t\t%d\n",        shell->mesh->Np*shell->mesh->Nelements);
    printf("Setup Time\t%g\n",      shell->times.setup.total);
    printf("Solve Time\t%g\n",      shell->times.solve);
    printf("Total Time\t%g\n",      shell->times.total);
    printf("Nodes/s\t\t%g\n", (1.0*shell->mesh->Np*shell->mesh->Nelements)/shell->times.solve);
    printf("Precon. memory\t%lld\n", shell->preconBytes);
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    printf("N\t\t%d\n",             shell->elliptic->mesh->N);
    printf("Np\t\t%d\n",            shell->elliptic->mesh->Np);
    printf("Nelements\t%d\n",       shell->elliptic->mesh->Nelements);
    printf("Ntotal\t\t%d\n",        shell->elliptic->mesh->Np*shell->elliptic->mesh->Nelements);
    printf("Setup Time\t%g\n",      shell->times.setup.total);
    printf("Solve Time\t%g\n",      shell->times.solve);
    printf("Total Time\t%g\n",      shell->times.total);
    printf("Nodes/s\t\t%g\n", (1.0*shell->elliptic->mesh->Np*shell->elliptic->mesh->Nelements)/shell->times.solve);
    printf("Precon. memory\t%lld\n", shell->preconBytes);
  }
  printf("------------------------\n");
  printf("H1 error:  %g\n", errH1);
  printf("L2 error:  %g\n", errL2);
 
  // close down MPI
  MPI_Finalize();

  return 0;
}
