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
#include "asbf.h"

int main(int argc, char **argv){
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

  // build asbf 
  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  asbf_t *asbf = asbfSetup(mesh, lambda, kernelInfo, options);
  
  // solve for asbf->q3D
  /*
  asbfSolve(asbf, options);

  // plot solution and compute error 
  if(asbf->elementType==QUADRILATERALS){
    dfloat errH1, errL2;
    // char fname[] = "sol";
    // asbfPlotVTU3D(asbf, fname, 0);
    asbfErrorHex3D(asbf, asbf->q3D, &errH1, &errL2);
    printf("%g, %g %%%% abs H1 error norm, L2 error norm\n", errH1, errL2);
  }

  exit(0);
  */

  /******************************/

  meshCheckHex3D(asbf->meshSEM);
  //meshPlotVTU3D(asbf->meshSEM, "semMesh", 0);

  occa::properties kernelInfoHex;
  kernelInfoHex["defines"].asObject();
  kernelInfoHex["includes"].asArray();
  kernelInfoHex["header"].asArray();
  kernelInfoHex["flags"].asObject();

  meshOccaSetup3D(asbf->meshSEM, options, kernelInfoHex);

  elliptic_t *hexSolver = new elliptic_t;
  asbf->meshSEM->Nfields = 1;
  
  hexSolver->mesh = asbf->meshSEM;
  hexSolver->options = asbf->options;
  hexSolver->dim = 3;
  hexSolver->elementType = HEXAHEDRA;

  int hexBCType[3] = {0, 1, 2};
  hexSolver->BCType = (int*)calloc(3, sizeof(int));
  memcpy(hexSolver->BCType, hexBCType, 3*sizeof(int));

  ellipticSolveSetup(hexSolver, asbf->lambda, kernelInfoHex);

  //---------------------------------------------------------------------------
  // Check symmetry by comparing y^TAx and x^TAy for random x, y.

  /*
  int Ntotal = hexSolver->mesh->Nelements*hexSolver->mesh->Np;
  dfloat *x  = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *y  = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *Ax = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *Ay = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  occa::memory o_v  = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));
  occa::memory o_Av = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));

  hexSolver->tau = 2*pow(hexSolver->mesh->N+1,2);

  srand48(67714070);

  for (int i = 0; i < Ntotal; i++) {
    x[i] = drand48();
    y[i] = drand48();
  }

  o_v.copyFrom(x);

  // feed o_v = S*G*o_v
  if(options.compareArgs("DISCRETIZATION", "CONTINUOUS")){
    ogsGatherScatterStart(o_v, ogsDfloat, ogsAdd, hexSolver->ogs);
    ogsGatherScatterFinish(o_v, ogsDfloat, ogsAdd, hexSolver->ogs);
  }

  ellipticOperator(hexSolver, asbf->lambda, o_v, o_Av, "double");
  o_Av.copyTo(Ax);
  double ytAx = 0.0;
  for (int i = 0; i < Ntotal; i++)
    ytAx += y[i]*Ax[i];

  o_v.copyFrom(y);

  // feed o_v = S*G*o_v
  if(options.compareArgs("DISCRETIZATION", "CONTINUOUS")){
    ogsGatherScatterStart(o_v, ogsDfloat, ogsAdd, hexSolver->ogs);
    ogsGatherScatterFinish(o_v, ogsDfloat, ogsAdd, hexSolver->ogs);
  }
  
  ellipticOperator(hexSolver, asbf->lambda, o_v, o_Av, "double");
  o_Av.copyTo(Ay);
  double xtAy = 0.0;
  for (int i = 0; i < Ntotal; i++)
    xtAy += x[i]*Ay[i];

  printf("ytAx = %.15f\n", ytAx);
  printf("xtAy = %.15f\n", xtAy);
  */

  //---------------------------------------------------------------------------

  hexSolver->x = (dfloat*)calloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np, sizeof(dfloat));
  hexSolver->r = (dfloat*)calloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np, sizeof(dfloat));

  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    for (int n = 0; n < hexSolver->mesh->Np; n++) {
      dfloat J = hexSolver->mesh->vgeo[hexSolver->mesh->Np*(e*hexSolver->mesh->Nvgeo + JID) + n];
      
      dlong ind = e*hexSolver->mesh->Np + n;
      dfloat xn = hexSolver->mesh->x[ind];
      dfloat yn = hexSolver->mesh->y[ind];
      dfloat zn = hexSolver->mesh->z[ind];

      dfloat r = sqrt(xn*xn + yn*yn + zn*zn);

      dfloat q, d2qdx2, d2qdy2, d2qdz2;
      dfloat r2mR2      = pow(r, 2.0) - pow(asbf->R, 2.0);
      dfloat r2m1       = pow(r, 2.0) - 1.0;
      dfloat twor2mR2m1 = 2.0*pow(r, 2.0) - pow(asbf->R, 2.0) - 1.0;
      
      q = sin(xn)*cos(yn)*exp(zn)*r2mR2*r2m1;
      d2qdx2 = cos(yn)*exp(zn)*((2.0*sin(xn) + 4.0*xn*cos(xn))*twor2mR2m1 + 8.0*xn*xn*sin(xn) - sin(xn)*r2mR2*r2m1);
      d2qdy2 = sin(xn)*exp(zn)*((2.0*cos(yn) - 4.0*yn*sin(yn))*twor2mR2m1 + 8.0*yn*yn*cos(yn) - cos(yn)*r2mR2*r2m1);
      d2qdz2 = sin(xn)*cos(yn)*exp(zn)*((2.0 + 4.0*zn)*twor2mR2m1 + 8.0*zn*zn + r2mR2*r2m1);

      hexSolver->r[ind] = J*(-d2qdx2 - d2qdy2 - d2qdz2 + asbf->lambda*q);
    }
  }
  
  if (options.compareArgs("BASIS","NODAL"))
    meshApplyElementMatrix(hexSolver->mesh,hexSolver->mesh->MM,hexSolver->r,hexSolver->r);

  hexSolver->o_r =
    hexSolver->mesh->device.malloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np*sizeof(dfloat), hexSolver->r);
  hexSolver->o_x =
    hexSolver->mesh->device.malloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np*sizeof(dfloat), hexSolver->x);

  if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ogsGatherScatter(hexSolver->o_r, ogsDfloat, ogsAdd, hexSolver->mesh->ogs);
    if (hexSolver->Nmasked)
      hexSolver->mesh->maskKernel(hexSolver->Nmasked, hexSolver->o_maskIds, hexSolver->o_r);
  }
  
  ellipticSolve(hexSolver, asbf->lambda, asbf->pTOL, hexSolver->o_r, hexSolver->o_x);

  hexSolver->o_x.copyTo(hexSolver->mesh->q);

  dfloat *gradx = (dfloat*)calloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np, 4*sizeof(dfloat));
  occa::memory o_gradx = hexSolver->mesh->device.malloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np*4*sizeof(dfloat));
  hexSolver->gradientKernel(hexSolver->mesh->Nelements, hexSolver->mesh->o_vgeo, hexSolver->mesh->o_D, hexSolver->o_x, o_gradx);
  o_gradx.copyTo(gradx);

  dfloat err = 0.0;
  dfloat graderrx = 0.0;
  dfloat graderry = 0.0;
  dfloat graderrz = 0.0;
  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    for (int i = 0; i < hexSolver->mesh->Np; i++) {
      dfloat x = hexSolver->mesh->x[e*hexSolver->mesh->Np + i];
      dfloat y = hexSolver->mesh->y[e*hexSolver->mesh->Np + i];
      dfloat z = hexSolver->mesh->z[e*hexSolver->mesh->Np + i];

      dfloat r          = sqrt(x*x + y*y + z*z);
      dfloat r2mR2      = pow(r, 2.0) - pow(asbf->R, 2.0);
      dfloat r2m1       = pow(r, 2.0) - 1.0;
      dfloat twor2mR2m1 = 2.0*pow(r, 2.0) - pow(asbf->R, 2.0) - 1.0;

      dfloat q = sin(x)*cos(y)*exp(z)*r2mR2*r2m1;
      dfloat dqdx = cos(y)*exp(z)*(2.0*x*sin(x)*twor2mR2m1 + cos(x)*r2mR2*r2m1); 
      dfloat dqdy = sin(x)*exp(z)*(2.0*y*cos(y)*twor2mR2m1 - sin(y)*r2mR2*r2m1);
      dfloat dqdz = sin(x)*cos(y)*exp(z)*(2.0*z*twor2mR2m1 + r2mR2*r2m1);

      dfloat thiserr = fabs(hexSolver->mesh->q[e*hexSolver->mesh->Np + i] - q);
      dfloat thisgraderrx = fabs(gradx[e*hexSolver->mesh->Np*4 + i*4 + 0] - dqdx);
      dfloat thisgraderry = fabs(gradx[e*hexSolver->mesh->Np*4 + i*4 + 1] - dqdy);
      dfloat thisgraderrz = fabs(gradx[e*hexSolver->mesh->Np*4 + i*4 + 2] - dqdz);

      //printf("error at (%g, %g, %g):  %g\n", x, y, z, thiserr);

      err = mymax(err, thiserr);
      graderrx = mymax(graderrx, thisgraderrx);
      graderry = mymax(graderry, thisgraderry);
      graderrz = mymax(graderrz, thisgraderrz);
    }
  }

  printf("error on grid:  %g\n", err);
  printf("x-gradient error on grid:  %g\n", graderrx);
  printf("y-gradient error on grid:  %g\n", graderry);
  printf("z-gradient error on grid:  %g\n", graderrz);

  //ellipticPlotVTUHex3D(hexSolver->mesh, "hexsol", 0);

  /******************************/
 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
