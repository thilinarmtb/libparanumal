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

  meshCheckHex3D(asbf->meshSEM);

  meshPlotVTU3D(asbf->meshSEM, "semMesh", 0);
  
  /*
  // solve for asbf->q3D
  asbfSolve(asbf, options);

  // plot solution and compute error 
  if(asbf->elementType==QUADRILATERALS){
    dfloat errH1, errL2;
    char fname[] = "sol";
    asbfPlotVTU3D(asbf, fname, 0);
    asbfErrorHex3D(asbf, asbf->q3D, &errH1, &errL2);
    printf("%g, %g %%%% abs H1 error norm, L2 error norm\n", errH1, errL2);
  }
  */

 /******************************/

  // maps nothing->nothing, Dirichlet->Dirichlet...
  int hexBCType[3] = {0,1,2};
  
  occa::properties kernelInfoHex;
  kernelInfoHex["defines"].asObject();
  kernelInfoHex["includes"].asArray();
  kernelInfoHex["header"].asArray();
  kernelInfoHex["flags"].asObject();

  meshOccaSetup3D(asbf->meshSEM, options, kernelInfoHex);

  elliptic_t *hexSolver = new elliptic_t[1];
  //(elliptic_t*)calloc(1, sizeof(elliptic_t));
  asbf->meshSEM->Nfields = 1;
  
  hexSolver->mesh = asbf->meshSEM;
  hexSolver->options = asbf->options;
  hexSolver->dim = 3;
  hexSolver->elementType = HEXAHEDRA;

  hexSolver->BCType = (int*) calloc(3,sizeof(int));
  memcpy(hexSolver->BCType,hexBCType,3*sizeof(int));

  ellipticSolveSetup(hexSolver, asbf->lambda, kernelInfoHex);

#if 0
  printf("EToB\n");
  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    printf("Element %d:  ", e);
    for (int f = 0; f < hexSolver->mesh->Nfaces; f++)
      printf("%d  ", hexSolver->mesh->EToB[e*hexSolver->mesh->Nfaces + f]);
    printf("\n");
  }

  printf("\n");

  printf("EToV\n");
  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    printf("Element %d:  ", e);
    for (int v = 0; v < hexSolver->mesh->Nverts; v++)
      printf("%d ", hexSolver->mesh->EToV[e*hexSolver->mesh->Nverts + v]);
    printf("\n");
  }

  printf("\n");

  printf("EToE\n");
  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    printf("Element %d:  ", e);
    for (int f = 0; f < hexSolver->mesh->Nfaces; f++)
      printf("%d  ", hexSolver->mesh->EToE[e*hexSolver->mesh->Nfaces + f]);
    printf("\n");
  }

  printf("\n");

  printf("EToF\n");
  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    printf("Element %d:  ", e);
    for (int f = 0; f < hexSolver->mesh->Nfaces; f++)
      printf("%d  ", hexSolver->mesh->EToF[e*hexSolver->mesh->Nfaces + f]);
    printf("\n");
  }

  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    for (int f = 0; f < hexSolver->mesh->Nfaces; f++) {
      for (int i = 0; i < hexSolver->mesh->Nfp; i++) {
        dlong base = hexSolver->mesh->Nsgeo*(hexSolver->mesh->Nfaces*hexSolver->mesh->Nfp*e + hexSolver->mesh->Nfp*f + i);
        if (((e == 297) && (f == 2)) || ((e == 346) && (f == 4))) {
          printf("Element %d, Face %d, Node %02d:  normal = (% 20.15f, % 20.15f, % 20.15f), surface Jacobian = % 20.15f\n",
                 e, f, i, hexSolver->mesh->sgeo[base + NXID], hexSolver->mesh->sgeo[base + NYID], hexSolver->mesh->sgeo[base + NZID],
                 hexSolver->mesh->sgeo[base + SJID]);
        }
      }
    }
  }

  exit(0);
#endif
  
  //---------------------------------------------------------------------------

  // Compute and dump the stiffness matrix to disk (this is too expensive).
  /*
  int Ntotal = hexSolver->mesh->Nelements*hexSolver->mesh->Np;
  dfloat *x  = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *Ax = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  occa::memory o_x  = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));
  occa::memory o_Ax = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));

  FILE *fout = fopen("stiffmat.dat", "w");

  int NN = 8;
  printf("Ntotal = %d\n", Ntotal);
  for (int i = 0; i < NN; i++) {
    if (i % 10000 == 0) {
      printf("i = %d\n", i);
    }

    x[i] = 1.0;
    if (i > 0)
      x[i - 1] = 0.0;

    o_x.copyFrom(x);
    ellipticOperator(hexSolver, asbf->lambda, o_x, o_Ax, "double");
    o_Ax.copyTo(Ax);

    for (int j = 0; j < NN; j++)
      printf("%.15f\n", Ax[j]);
    printf("\n");

    //for (int j = 0; j < Ntotal; j++) {
    //  if (fabs(Ax[j]) > 1.0e-14)
    //    fprintf(fout, "%d %d %.16f\n", j, i, Ax[j]);
    //}
  }

  free(x);
  free(Ax);
  fclose(fout);
  */

  // Check symmetry by comparing y^TAx and x^TAy for random x, y.
  int Ntotal = hexSolver->mesh->Nelements*hexSolver->mesh->Np;
  dfloat *x  = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *y  = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *Ax = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  dfloat *Ay = (dfloat*)calloc(Ntotal, sizeof(dfloat));
  occa::memory o_v  = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));
  occa::memory o_Av = hexSolver->mesh->device.malloc(Ntotal*sizeof(dfloat));

  hexSolver->tau = 1000.0;

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

  //---------------------------------------------------------------------------

  /*
  hlong Nall = hexSolver->mesh->Np*(hexSolver->mesh->Nelements+hexSolver->mesh->totalHaloPairs);
  dfloat *foo = (dfloat*) malloc(Nall*sizeof(dfloat));
  hexSolver->o_invDegree.copyTo(foo);
  */

  hexSolver->x = (dfloat*)calloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np, sizeof(dfloat));
  hexSolver->r = (dfloat*)calloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np, sizeof(dfloat));

  for (int e = 0; e < hexSolver->mesh->Nelements; e++) {
    for (int n = 0; n < hexSolver->mesh->Np; n++) {
      dfloat J = hexSolver->mesh->vgeo[hexSolver->mesh->Np*(e*hexSolver->mesh->Nvgeo + JID) + n];

      dlong ind = e*hexSolver->mesh->Np + n;
      dfloat xn = hexSolver->mesh->x[ind];
      dfloat yn = hexSolver->mesh->y[ind];
      dfloat zn = hexSolver->mesh->z[ind];
      dfloat r = sqrt(zn*zn + yn*yn + zn*zn);

      dfloat q, d2qdx2, d2qdy2, d2qdz2;
      dfloat r2mR2      = pow(r, 2.0) - pow(asbf->R, 2.0);
      dfloat r2m1       = pow(r, 2.0) - 1.0;
      dfloat twor2mR2m1 = 2.0*pow(r, 2.0) - pow(asbf->R, 2.0) - 1.0;
      
      q = sin(xn)*cos(yn)*exp(zn)*r2mR2*r2m1;
      d2qdx2 = cos(yn)*exp(zn)*((2.0*sin(xn) + 4.0*xn*cos(xn))*twor2mR2m1 + 8.0*xn*xn*sin(xn) - sin(xn)*r2mR2*r2m1);
      d2qdy2 = sin(xn)*exp(zn)*((2.0*cos(yn) - 4.0*yn*sin(yn))*twor2mR2m1 + 8.0*yn*yn*cos(yn) - cos(yn)*r2mR2*r2m1);
      d2qdz2 = sin(xn)*cos(yn)*exp(zn)*((2.0 + 4.0*zn)*twor2mR2m1 + 8.0*zn*zn + r2mR2*r2m1);
      hexSolver->r[ind] = -d2qdx2 - d2qdy2 - d2qdz2 + asbf->lambda*q;
    }
  }
  
  hexSolver->o_r = hexSolver->mesh->device.malloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np*sizeof(dfloat), hexSolver->r);
  hexSolver->o_x = hexSolver->mesh->device.malloc(hexSolver->mesh->Nelements*hexSolver->mesh->Np*sizeof(dfloat), hexSolver->x);

  ellipticSolve(hexSolver, asbf->lambda, asbf->pTOL, hexSolver->o_r, hexSolver->o_x);

  hexSolver->o_x.copyTo(hexSolver->mesh->q);

 /******************************/
 
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
