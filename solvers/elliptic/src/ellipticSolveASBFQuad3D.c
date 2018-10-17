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

#include "elliptic.h"

int ellipticSolveASBFQuad3D(elliptic_t *elliptic,
			    dfloat lambda,
			    dfloat tol,
			    occa::memory &o_r,
			    occa::memory &o_q){

  mesh_t *mesh = elliptic->mesh;
  
  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/quadrilateralN%02d.dat", mesh->N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  int Nrows, Ncols;
  
  readDfloatArray(fp, "ASBF EIGENVALUES",
		  &(elliptic->asbfEigenvalues),&(elliptic->asbfNmodes), &(Ncols));
  
  readDfloatArray(fp, "ASBF QUADRATURE VANDERMONDE",
		  &(elliptic->asbfBquad),&(elliptic->asbfNquad), &(elliptic->asbfNmodes));

  readDfloatArray(fp, "ASBF QUADRATURE NODES",
		  &(elliptic->asbfRquad),&(elliptic->asbfNquad), &Ncols);

  readDfloatArray(fp, "ASBF QUADRATURE WEIGHTS",
		  &(elliptic->asbfWquad),&(elliptic->asbfNquad), &Ncols);
  
  readDfloatArray(fp, "ASBF GLL VANDERMONDE",
		  &(elliptic->asbfBgll),&(elliptic->asbfNgll), &(elliptic->asbfNmodes));
  
  readDfloatArray(fp, "ASBF GLL NODES",
		  &(elliptic->asbfRgll),&(elliptic->asbfNgll), &Ncols);

  dlong Ntotal = mesh->Np*mesh->Nelements;
  dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
  dlong Nall   = Ntotal + Nhalo;
  
  dfloat *r3D = (dfloat*) calloc(elliptic->asbfNmodes*Nall, sizeof(dfloat));
  dfloat *q3D = (dfloat*) calloc(elliptic->asbfNmodes*Nall, sizeof(dfloat));
  dfloat *f   = (dfloat*) calloc(elliptic->asbfNquad, sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat xbase = mesh->x[e*mesh->Np+n];
      dfloat ybase = mesh->y[e*mesh->Np+n];
      dfloat zbase = mesh->z[e*mesh->Np+n];

      dfloat JW = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JWID) + n];
      
      for(int g=0;g<elliptic->asbfNquad;++g){
	
	dfloat Rg = elliptic->asbfRquad[g];

	// stretch coordinates 
	dfloat xg = Rg*xbase;
	dfloat yg = Rg*ybase;
	dfloat zg = Rg*zbase;

	// evaluate rhs at asbf quadrature for each surface node
	f[g] = sin(M_PI*xg)*sin(M_PI*yg)*sin(M_PI*zg);
      }

      // integrate f against asbf modes
      for(int m=0;m<elliptic->asbfNmodes;++m){
	dfloat fhatm = 0;
	for(int i=0;i<elliptic->asbfNquad;++i){
	  fhatm += elliptic->asbfBquad[m + i*elliptic->asbfNmodes]*elliptic->asbfWquad[i]*f[i];
	}

	// scale by surface weight
	r3D[e*mesh->Np + n + m*Nall] = JW*fhatm;
      }
    }
  }

  for(int m=0;m<elliptic->asbfNmodes;++m){

    o_r.copyFrom(r3D + Nall*m);
    
    dfloat lambdam = lambda + elliptic->asbfEigenvalues[m];

    ellipticSolve(elliptic, lambdam, tol, o_r, o_q);
    
    o_q.copyTo(q3D + Nall*m);
  }

  mesh_t *meshSEM = (mesh_t*) calloc(1, sizeof(mesh_t*));
  
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
  
  meshSEM->x = (dfloat*) calloc(mesh->Nq*Nall, sizeof(dfloat));
  meshSEM->y = (dfloat*) calloc(mesh->Nq*Nall, sizeof(dfloat));
  meshSEM->z = (dfloat*) calloc(mesh->Nq*Nall, sizeof(dfloat));
  meshSEM->q = (dfloat*) calloc(mesh->Nq*Nall, sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat xbase = mesh->x[e*mesh->Np+n];
      dfloat ybase = mesh->y[e*mesh->Np+n];
      dfloat zbase = mesh->z[e*mesh->Np+n];

      for(int g=0;g<mesh->Nq;++g){
	dfloat Rg = elliptic->asbfRgll[g];
	
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
	for(int i=0;i<elliptic->asbfNmodes;++i){
	  qg += elliptic->asbfBgll[i + g*elliptic->asbfNmodes]*q3D[e*mesh->Np+n+i*Nall];
	}
	// assume Nfields=1
	meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
      }
    }
  }
  
  ellipticPlotVTUHex3D(meshSEM, "foo", 0);
}

