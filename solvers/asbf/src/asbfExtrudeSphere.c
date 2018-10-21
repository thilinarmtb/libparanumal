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

void asbfExtrudeSphere(asbf_t *asbf){

  mesh_t *mesh = asbf->mesh;
  
  mesh_t *meshSEM = (mesh_t*) calloc(1, sizeof(mesh_t));
  
  meshSEM->N = mesh->N;
  meshSEM->Nelements = mesh->Nelements;
  meshSEM->dim = 3;
  meshSEM->Nverts = 8;
  meshSEM->Nfaces = 6;
  meshSEM->NfaceVertices = 4;
  
  // yuck --->
  int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};
  meshSEM->faceVertices =
    (int*) calloc(meshSEM->NfaceVertices*meshSEM->Nfaces, sizeof(int));
  memcpy(meshSEM->faceVertices, faceVertices[0], meshSEM->NfaceVertices*meshSEM->Nfaces*sizeof(int));
  // <----
  
  meshSEM->rank = mesh->rank;
  meshSEM->size = mesh->size;
  meshSEM->comm = mesh->comm;
  
  // clone and shift nodes
  meshSEM->Nnodes = 2*mesh->Nnodes;
  meshSEM->EToV = (hlong*)  calloc(meshSEM->Nelements*meshSEM->Nverts, sizeof(hlong));
  meshSEM->EX   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));
  meshSEM->EY   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));
  meshSEM->EZ   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));
  
  for(dlong e=0;e<meshSEM->Nelements;++e){
    for(int n=0;n<mesh->Nverts;++n){
      hlong id   = e*mesh->Nverts+n;
      hlong vid1 = mesh->EToV[id];
      hlong vid2 = vid1 + mesh->Nnodes;

      dfloat Rg1 = asbf->asbfRgll[0];
      hlong id1  = e*meshSEM->Nverts+n;
      meshSEM->EToV[id1] = vid1;
      meshSEM->EX[id1] = Rg1*mesh->EX[id];
      meshSEM->EY[id1] = Rg1*mesh->EY[id];
      meshSEM->EZ[id1] = Rg1*mesh->EZ[id];
      
      dfloat Rg2 = asbf->asbfRgll[mesh->Nq-1];
      hlong id2  = e*meshSEM->Nverts+n+mesh->Nverts;
      meshSEM->EToV[id2] = vid2;
      meshSEM->EX[id2] = Rg2*mesh->EX[id];
      meshSEM->EY[id2] = Rg2*mesh->EY[id];
      meshSEM->EZ[id2] = Rg2*mesh->EZ[id];
    }
  }
  
  meshSEM->NboundaryFaces = 2*mesh->Nelements;
  meshSEM->boundaryInfo = (hlong*) calloc(meshSEM->NboundaryFaces*(meshSEM->NfaceVertices+1), sizeof(hlong));
  hlong bcnt = 0;
  for(hlong e=0;e<meshSEM->Nelements;++e){
    // inner sphere surface
    meshSEM->boundaryInfo[bcnt*5+0] = 1; // DIRICHLET 
    meshSEM->boundaryInfo[bcnt*5+1] = meshSEM->EToV[e*meshSEM->Nverts+0];
    meshSEM->boundaryInfo[bcnt*5+2] = meshSEM->EToV[e*meshSEM->Nverts+1];
    meshSEM->boundaryInfo[bcnt*5+3] = meshSEM->EToV[e*meshSEM->Nverts+2];
    meshSEM->boundaryInfo[bcnt*5+4] = meshSEM->EToV[e*meshSEM->Nverts+3];
    ++bcnt;
    
    // outer sphere surface 
    meshSEM->boundaryInfo[bcnt*5+0] = 1; // DIRICHLET ?
    meshSEM->boundaryInfo[bcnt*5+1] = meshSEM->EToV[e*meshSEM->Nverts+0+mesh->Nverts];
    meshSEM->boundaryInfo[bcnt*5+2] = meshSEM->EToV[e*meshSEM->Nverts+1+mesh->Nverts];
    meshSEM->boundaryInfo[bcnt*5+3] = meshSEM->EToV[e*meshSEM->Nverts+2+mesh->Nverts];
    meshSEM->boundaryInfo[bcnt*5+4] = meshSEM->EToV[e*meshSEM->Nverts+3+mesh->Nverts];
    ++bcnt;
  }
  
  // connect elements using parallel sort
  meshParallelConnect(meshSEM);
  
  // print out connectivity statistics
  meshPartitionStatistics(meshSEM);
  
  // connect elements to boundary faces
  meshConnectBoundary(meshSEM);
  
  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(meshSEM, meshSEM->N);
  
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
    }
  }
  
  // compute geometric factors
  meshGeometricFactorsHex3D(meshSEM);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(meshSEM);
  
  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(meshSEM);
  
  // compute surface geofacs (including halo)
  meshSurfaceGeometricFactorsHex3D(meshSEM);
  
  // global nodes
  meshParallelConnectNodes(meshSEM);
  
  asbf->meshSEM = meshSEM;
}
