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
#include "omp.h"
#include <unistd.h>

asbf_t *asbfSetup(mesh_t *mesh, setupAide options){

  asbf_t *asbf = (asbf_t*) calloc(1, sizeof(asbf_t));
  asbf->mesh = mesh;
  asbf->options = options;

  options.getArgs("MESH DIMENSION", asbf->dim);
  options.getArgs("ELEMENT TYPE", asbf->elementType);

  mesh->Nfields = 1;

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/solvers/asbf/data/asbfN%02d.dat", mesh->N);

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  int Nrows, Ncols;

  readDfloatArray(fp, "ASBF EIGENVALUES",
		  &(asbf->asbfEigenvalues),&(asbf->asbfNmodes), &(Ncols));
  readDfloatArray(fp, "ASBF QUADRATURE VANDERMONDE",
		  &(asbf->asbfBquad),&(asbf->asbfNquad), &(asbf->asbfNmodes));
  readDfloatArray(fp, "ASBF QUADRATURE NODES",
		  &(asbf->asbfRquad),&(asbf->asbfNquad), &Ncols);
  readDfloatArray(fp, "ASBF QUADRATURE WEIGHTS",
		  &(asbf->asbfWquad),&(asbf->asbfNquad), &Ncols);
  readDfloatArray(fp, "ASBF GLL VANDERMONDE",
		  &(asbf->asbfBgll),&(asbf->asbfNgll), &(asbf->asbfNmodes));
  readDfloatArray(fp, "ASBF GLL NODES",
		  &(asbf->asbfRgll),&(asbf->asbfNgll), &Ncols);

  fclose(fp);

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
  asbf->Ntotal = Nlocal + Nhalo;

  asbf->r3D = (dfloat*) calloc(asbf->asbfNmodes*asbf->Ntotal, sizeof(dfloat));
  asbf->q3D = (dfloat*) calloc(asbf->asbfNmodes*asbf->Ntotal, sizeof(dfloat));
  asbf->f   = (dfloat*) calloc(asbf->asbfNquad, sizeof(dfloat));
  asbf->r   = (dfloat*) calloc(asbf->Ntotal, sizeof(dfloat));
  asbf->x   = (dfloat*) calloc(asbf->Ntotal, sizeof(dfloat));

  for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){

        dfloat xbase = mesh->x[e*mesh->Np+n];
        dfloat ybase = mesh->y[e*mesh->Np+n];
        dfloat zbase = mesh->z[e*mesh->Np+n];

        dfloat JW = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JWID) + n];

        for(int g=0;g<asbf->asbfNquad;++g){

          dfloat Rg = asbf->asbfRquad[g];

          // stretch coordinates
          dfloat xg = Rg*xbase;
          dfloat yg = Rg*ybase;
          dfloat zg = Rg*zbase;

          // evaluate rhs at asbf quadrature for each surface node
          asbf->f[g] = sin(M_PI*xg)*sin(M_PI*yg)*sin(M_PI*zg);
        }

        // integrate f against asbf modes
        for(int m=0;m<asbf->asbfNmodes;++m){
          dfloat fhatm = 0;
          for(int i=0;i<asbf->asbfNquad;++i){
            fhatm += asbf->asbfBquad[m + i*asbf->asbfNmodes]*asbf->asbfWquad[i]*asbf->f[i];
          }

          // scale by surface weight
          asbf->r3D[e*mesh->Np + n + m*asbf->Ntotal] = JW*fhatm;
        }
      }
    }

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

   if(asbf->dim==3){
     if(asbf->elementType != QUADRILATERALS)
       meshOccaSetup3D(mesh, options, kernelInfo);
     else
       meshOccaSetupQuad3D(mesh, options, kernelInfo);
   }
   else
     meshOccaSetup2D(mesh, options, kernelInfo);

  occa::properties kernelInfoP  = kernelInfo;

  asbf->o_r   = mesh->device.malloc(asbf->Ntotal*sizeof(dfloat), asbf->r);
  asbf->o_x   = mesh->device.malloc(asbf->Ntotal*sizeof(dfloat), asbf->x);

  asbf->pOptions = options;
  asbf->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("KRYLOV SOLVER"));
  asbf->pOptions.setArgs("DISCRETIZATION",       options.getArgs("DISCRETIZATION"));
  asbf->pOptions.setArgs("BASIS",                options.getArgs("BASIS"));
  asbf->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRECONDITIONER"));
  asbf->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("MULTIGRID COARSENING"));
  asbf->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("MULTIGRID SMOOTHER"));
  asbf->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PARALMOND CYCLE"));
  asbf->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PARALMOND SMOOTHER"));
  asbf->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PARALMOND PARTITION"));

  options.getArgs("LAMBDA", asbf->lambda);

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int pBCType[7] = {0,1,1,2,1,1,1}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances
  asbf->pTOL = 1E-8;

  asbf->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  asbf->pSolver->mesh = mesh;
  asbf->pSolver->options = asbf->pOptions;
  asbf->pSolver->dim = asbf->dim;
  asbf->pSolver->elementType = asbf->elementType;
  asbf->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(asbf->pSolver->BCType,pBCType,7*sizeof(int));

  ellipticSolveSetup(asbf->pSolver, asbf->lambda, kernelInfoP);

  mesh_t *meshSEM = (mesh_t*) calloc(1, sizeof(mesh_t));

  meshSEM->Nelements = mesh->Nelements;
  meshSEM->dim = 3;
  meshSEM->Nverts = 8;
  meshSEM->Nfaces = 6;
  meshSEM->NfaceVertices = 4;

  // yuck --->
  int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(int));
  // <----

  meshSEM->rank = mesh->rank;
  meshSEM->size = mesh->size;
  meshSEM->comm = mesh->comm;

  meshLoadReferenceNodesHex3D(meshSEM, mesh->N);

  meshSEM->x = (dfloat*) calloc(mesh->Nq*asbf->Ntotal, sizeof(dfloat));
  meshSEM->y = (dfloat*) calloc(mesh->Nq*asbf->Ntotal, sizeof(dfloat));
  meshSEM->z = (dfloat*) calloc(mesh->Nq*asbf->Ntotal, sizeof(dfloat));
  meshSEM->q = (dfloat*) calloc(mesh->Nq*asbf->Ntotal, sizeof(dfloat));

  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat xbase = mesh->x[e*mesh->Np+n];
      dfloat ybase = mesh->y[e*mesh->Np+n];
      dfloat zbase = mesh->z[e*mesh->Np+n];

      for(int g=0;g<mesh->Nq;++g){
	dfloat Rg = asbf->asbfRgll[g];

	// stretch coordinates
	dfloat xg = Rg*xbase;
	dfloat yg = Rg*ybase;
	dfloat zg = Rg*zbase;
	meshSEM->x[e*meshSEM->Np+g*mesh->Np+n] = xg;
	meshSEM->y[e*meshSEM->Np+g*mesh->Np+n] = yg;
	meshSEM->z[e*meshSEM->Np+g*mesh->Np+n] = zg;
      }

      // interpolate from asbf to gll nodes
      for(int g=0;g<mesh->Nq;++g){
	dfloat qg = 0;
	for(int i=0;i<asbf->asbfNmodes;++i){
	  qg += asbf->asbfBgll[i + g*asbf->asbfNmodes]*asbf->q3D[(e*mesh->Np+n)+i*asbf->Ntotal];
	}
	// assume Nfields=1
	meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
      }
    }
  }

  asbf->meshSEM = meshSEM;

  return asbf;
}
