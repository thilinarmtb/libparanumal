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

#include "cuda.h"
#include "cudaProfiler.h"
#include "elliptic.h"
#include "timer.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  // set up mesh stuff

  int dim, elementType;

  string fileName = argv[2];
  int    N        = atoi(argv[3]);
  int    maxiter  = atoi(argv[4]);
  int 	 ppn	  = atoi(argv[5]);

  int NTEST      = 1; // number of test for BP5 
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  // set up mesh
  mesh_t *mesh = meshSetup((char*) fileName.c_str(), N, options);

  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  // set up
  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  elliptic_t *elliptic = ellipticSetup(mesh, lambda, kernelInfo, options);

#if(TIMER) 
  timer *profiler = elliptic->profiler; 
  profiler->tic("PCG");
#endif

  // do a hot start
  pcgBP5(elliptic, lambda, elliptic->o_r, elliptic->o_x, 100);

#if 1
  double start = 0.0, end =0.0;
  mesh->device.finish();
  MPI_Barrier(mesh->comm);
  start = MPI_Wtime();
#else
  occa::streamTag startTag = mesh->device.tagStream();
#endif

  for(int tst=0; tst<NTEST; tst++){
    pcgBP5(elliptic, lambda, elliptic->o_r, elliptic->o_x, maxiter);
  }


#if 1
  mesh->device.finish();
  MPI_Barrier(mesh->comm);
  end = MPI_Wtime();
  double localElapsed = end-start;
#else
  occa::streamTag stopTag = mesh->device.tagStream();
  double localElapsed = mesh->device.timeBetween(startTag, stopTag);
#endif

  localElapsed /= NTEST; // Average time for each PCG solve; 

  int size = mesh->size;

  hlong   localDofs     = (hlong) mesh->Np*mesh->Nelements;
  hlong   localElements = (hlong) mesh->Nelements;
  double  globalElapsed;
  hlong   globalDofs;
  hlong   globalElements;

  MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE,  MPI_MAX, 0, mesh->comm );
  MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );
  MPI_Reduce(&localElements,&globalElements,1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );

  globalElapsed /= maxiter;
  if (mesh->rank==0){
    printf("libP,occa,BP5,%d,%d,%d,%lld,%g\n",
	   mesh->N,mesh->size,ppn,globalElements,globalElapsed);
  }

  // close down MPI
  MPI_Finalize();

  return 0;
}
