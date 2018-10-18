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

int asbfSolve(asbf_t *asbf, setupAide options)
{
  mesh_t *mesh         = asbf->mesh;
  elliptic_t *elliptic = asbf->pSolver;

  for(int m=0;m<asbf->asbfNmodes;++m){
    asbf->o_r.copyFrom(asbf->r3D + asbf->Ntotal*m);
    dfloat lambdam = asbf->lambda + asbf->asbfEigenvalues[m];
    ellipticSolve(elliptic, lambdam, asbf->pTOL, asbf->o_r, asbf->o_x);
    asbf->o_x.copyTo(asbf->q3D + asbf->Ntotal*m);
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
	  qg += asbf->asbfBgll[i + g*asbf->asbfNmodes]*asbf->q3D[e*mesh->Np+n+i*asbf->Ntotal];
	}
	// assume Nfields=1
	meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
      }
    }
  }


  char fname[] = "foo";
  ellipticPlotVTUHex3D(meshSEM, fname, 0);
}
