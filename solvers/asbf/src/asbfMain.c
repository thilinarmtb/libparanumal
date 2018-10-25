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

#include "asbf.h"

int main(int argc, char **argv){

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

  // build asbf 
  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  asbf_t *asbf = asbfSetup(mesh, lambda, kernelInfo, options);

  // solve for asbf->meshSEM->q
  asbfSolve(asbf, options);

  // compute error 
  if(asbf->elementType==QUADRILATERALS){

    char fname[] = "sol";
    asbfPlotVTU3D(asbf, fname, 0);
    //    ellipticPlotVTUHex3D(asbf->meshSEM, fname, 0);
    
    asbfErrorHex3D(asbf, asbf->q3D);
  }
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
