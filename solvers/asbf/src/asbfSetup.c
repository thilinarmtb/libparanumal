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

  options.getArgs("LAMBDA", asbf->lambda);
  
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

      //      dfloat JW = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JWID) + n];
      dfloat J;

      if(asbf->elementType==QUADRILATERALS)
	J = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JID) + n];
      else
	J = mesh->vgeo[e*mesh->Nvgeo + JID];
      
      for(int g=0;g<asbf->asbfNquad;++g){
	
	dfloat Rg = asbf->asbfRquad[g];
	
	// stretch coordinates
	dfloat xg = Rg*xbase;
	dfloat yg = Rg*ybase;
	dfloat zg = Rg*zbase;

	// evaluate rhs at asbf quadrature for each surface node
	//asbf->f[g] = sin(M_PI*xg)*sin(M_PI*yg)*sin(M_PI*zg);

	dfloat k1 = 6.283185307179586;
	dfloat k2 = 18.849555921538759;
	dfloat k3 = 25.132741228718345;
	dfloat r = sqrt(xg*xg + yg*yg + zg*zg);
	asbf->f[g] = (k1 + asbf->lambda/k1)*sin(k1*r)/r
	  + (k2 + asbf->lambda/k2)*sin(k2*r)/r
	  + (k3 + asbf->lambda/k3)*sin(k3*r)/r;

	//	asbf->f[g] = sin(M_PI*xg)*sin(M_PI*yg)*sin(M_PI*zg);
      }

      // integrate f against asbf modes
      for(int m=0;m<asbf->asbfNmodes;++m){
	dfloat fhatm = 0;
	for(int i=0;i<asbf->asbfNquad;++i){
	  fhatm += asbf->asbfBquad[m + i*asbf->asbfNmodes]*asbf->asbfWquad[i]*asbf->f[i];
	}

	// scale by surface weight
	//	asbf->r3D[e*mesh->Np + n + m*asbf->Ntotal] = JW*fhatm;
	asbf->r3D[e*mesh->Np + n + m*asbf->Ntotal] = J*fhatm;
      }
    }
  }
  
  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(asbf->dim==3){
    if(asbf->elementType == QUADRILATERALS){
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
    }else if(asbf->elementType == TRIANGLES){
      meshOccaSetupTri3D(mesh, options, kernelInfo);
    }else{
      meshOccaSetup3D(mesh, options, kernelInfo);
    }
  }
  else
    meshOccaSetup2D(mesh, options, kernelInfo);
  
  kernelInfo["defines/p_asbfNmodes"] = asbf->asbfNmodes;
  kernelInfo["defines/p_asbfNquad"] = asbf->asbfNquad;
  kernelInfo["defines/p_asbfNgll"] = asbf->asbfNgll;
   
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

  
  mesh_t *meshSEM = NULL;

  if(asbf->elementType==QUADRILATERALS){
    meshSEM = (mesh_t*) calloc(1, sizeof(mesh_t));

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
  }  

  asbf->meshSEM = meshSEM;

  // OKL kernels specific to asbf
  asbf->asbfReconstructKernel =
    mesh->device.buildKernel(DASBF "/okl/asbfReconstructHex3D.okl",
			     "asbfReconstructHex3D",
			     kernelInfo);

  return asbf;
}
