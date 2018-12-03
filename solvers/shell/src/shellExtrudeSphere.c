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

#include "shell.h"

void shellExtrudeSphere(shell_t *shell)
{
  mesh_t *mesh = shell->mesh;

  mesh_t *meshSEM = new mesh_t();

  meshSEM->N = mesh->N;
  meshSEM->Nelements = shell->Nradelements*mesh->Nelements;
  meshSEM->dim = 3;
  meshSEM->Nverts = 8;
  meshSEM->Nfaces = 6;
  meshSEM->NfaceVertices = 4;
  meshSEM->Nfields = 1;

  // yuck --->
  int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};
  meshSEM->faceVertices =
    (int*) calloc(meshSEM->NfaceVertices*meshSEM->Nfaces, sizeof(int));
  memcpy(meshSEM->faceVertices, faceVertices[0], meshSEM->NfaceVertices*meshSEM->Nfaces*sizeof(int));
  // <----

  meshSEM->rank = mesh->rank;
  meshSEM->size = mesh->size;
  MPI_Comm_dup(mesh->comm, &meshSEM->comm);

  int maxNode = 0;
  for(int n=0;n<mesh->Nelements*mesh->Nverts;++n){
    maxNode = mymax(maxNode, mesh->EToV[n]);
  }
  
  // clone and shift nodes
  printf("mesh->Nnodes = %d and maxNode = %d\n", mesh->Nnodes, maxNode);
  meshSEM->Nnodes = (shell->Nradelements + 1)*mesh->Nnodes;
  meshSEM->EToV = (hlong*)  calloc(meshSEM->Nelements*meshSEM->Nverts, sizeof(hlong));
  meshSEM->EX   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));
  meshSEM->EY   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));
  meshSEM->EZ   = (dfloat*) calloc(meshSEM->Nverts*meshSEM->Nelements, sizeof(dfloat));

  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int i = 0; i < shell->Nradelements; i++) {
      dlong eSEM = e*shell->Nradelements + i;
      dfloat rBot = shell->Rbreaks[i];
      dfloat rTop = shell->Rbreaks[i + 1];

      for (int n = 0; n < mesh->Nverts; n++) {
        hlong id = e*mesh->Nverts + n;
        hlong vid = mesh->EToV[id];

        hlong id1 = eSEM*meshSEM->Nverts + n;
        hlong vid1 = vid + i*mesh->Nnodes;
        meshSEM->EToV[id1] = vid1;
        meshSEM->EX[id1] = rBot*mesh->EX[id];
        meshSEM->EY[id1] = rBot*mesh->EY[id];
        meshSEM->EZ[id1] = rBot*mesh->EZ[id];

        hlong id2 = eSEM*meshSEM->Nverts + n + mesh->Nverts;
        hlong vid2 = vid + (i + 1)*mesh->Nnodes;
        meshSEM->EToV[id2] = vid2;
        meshSEM->EX[id2] = rTop*mesh->EX[id];
        meshSEM->EY[id2] = rTop*mesh->EY[id];
        meshSEM->EZ[id2] = rTop*mesh->EZ[id];
      }
    }
  }

  meshSEM->NboundaryFaces = 2*mesh->Nelements;
  meshSEM->boundaryInfo = (hlong*)calloc(meshSEM->NboundaryFaces*(1 + meshSEM->NfaceVertices), sizeof(hlong));
  hlong bcnt = 0;

  for (hlong e = 0; e < mesh->Nelements; e++) {
    // Inner sphere surface.
    hlong eSEMInner = e*shell->Nradelements;
    hlong idInner = eSEMInner*meshSEM->Nverts;
    meshSEM->boundaryInfo[bcnt*5 + 0] = shell->innerBC;
    meshSEM->boundaryInfo[bcnt*5 + 1] = meshSEM->EToV[idInner + 0];
    meshSEM->boundaryInfo[bcnt*5 + 2] = meshSEM->EToV[idInner + 1];
    meshSEM->boundaryInfo[bcnt*5 + 3] = meshSEM->EToV[idInner + 2];
    meshSEM->boundaryInfo[bcnt*5 + 4] = meshSEM->EToV[idInner + 3];
    bcnt++;

    // Outer sphere surface.
    hlong eSEMOuter = (e + 1)*shell->Nradelements - 1;
    hlong idOuter = eSEMOuter*meshSEM->Nverts + mesh->Nverts;
    meshSEM->boundaryInfo[bcnt*5 + 0] = shell->outerBC;
    meshSEM->boundaryInfo[bcnt*5 + 1] = meshSEM->EToV[idOuter + 0];
    meshSEM->boundaryInfo[bcnt*5 + 2] = meshSEM->EToV[idOuter + 1];
    meshSEM->boundaryInfo[bcnt*5 + 3] = meshSEM->EToV[idOuter + 2];
    meshSEM->boundaryInfo[bcnt*5 + 4] = meshSEM->EToV[idOuter + 3];
    bcnt++;
  }

  // connect elements using parallel sort
  meshParallelConnect(meshSEM);

  // print out connectivity statistics
  meshPartitionStatistics(meshSEM);

  // connect elements to boundary faces
  meshConnectBoundary(meshSEM);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(meshSEM, meshSEM->N);

  meshSEM->x = (dfloat*)calloc(mesh->Nq*shell->Ntotal*shell->Nradelements, sizeof(dfloat));
  meshSEM->y = (dfloat*)calloc(mesh->Nq*shell->Ntotal*shell->Nradelements, sizeof(dfloat));
  meshSEM->z = (dfloat*)calloc(mesh->Nq*shell->Ntotal*shell->Nradelements, sizeof(dfloat));
  meshSEM->q = (dfloat*)calloc(mesh->Nq*shell->Ntotal*shell->Nradelements, sizeof(dfloat));

  for (int e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) {
      dfloat xbase = mesh->x[e*mesh->Np + n];
      dfloat ybase = mesh->y[e*mesh->Np + n];
      dfloat zbase = mesh->z[e*mesh->Np + n];

      for (int i = 0; i < shell->Nradelements; i++) {
        int eSEM = e*shell->Nradelements + i;
        dfloat J = (shell->Rbreaks[i + 1] - shell->Rbreaks[i])/2.0;

        for (int j = 0; j < mesh->Nq; j++) {
          dfloat r = (shell->Rgll[j] + 1.0)*J + shell->Rbreaks[i];
          meshSEM->x[eSEM*meshSEM->Np + j*mesh->Np + n] = r*xbase;
          meshSEM->y[eSEM*meshSEM->Np + j*mesh->Np + n] = r*ybase;
          meshSEM->z[eSEM*meshSEM->Np + j*mesh->Np + n] = r*zbase;
        }
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

  shell->meshSEM = meshSEM;
}
