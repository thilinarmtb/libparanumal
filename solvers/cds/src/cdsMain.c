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

#include "cds.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./cdsMain setupfile\n");
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
  switch(elementType){
  case TRIANGLES:
    mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
  case QUADRILATERALS:{
    if(dim==2)
      mesh = meshSetupQuad2D((char*)fileName.c_str(), N);
  else{
    dfloat radius = 1;
    options.getArgs("SPHERE RADIUS", radius);
    mesh = meshSetupQuad3D((char*)fileName.c_str(), N, radius);
  }
  break;
  }
  case TETRAHEDRA:
    mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  cds_t *cds = cdsSetup(mesh,options);

   if (cds->options.compareArgs("TIME INTEGRATOR", "ARK")){
   printf("ARK is not implemented yet\n");
    exit(EXIT_FAILURE);
  }else if (cds->options.compareArgs("TIME INTEGRATOR", "EXTBDF")){
     if(mesh->rank==0) printf("Running EXTBDF Scheme\n");
     cdsRunEXTBDF(cds);
  }
  // close down MPI
  MPI_Finalize();
  exit(0);
  return 0;
}
