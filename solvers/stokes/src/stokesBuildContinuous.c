/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

#include "stokes.h"


// compare on global indices
int stokesParallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

void stokesBuildContinuousQuad2D(stokes_t *stokes, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);
void stokesBuildContinuousHex3D (stokes_t *stokes, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);



void stokesBuildContinuous(stokes_t *stokes, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

  switch(stokes->elementType){
  case QUADRILATERALS:{ // not 3d quads yet
    stokesBuildContinuousQuad2D(stokes, lambda, A, nnz, ogs, globalStarts);
    break;
  }
  default:
    printf("THIS ELEMENT NOT SUPPORTED BY stokesBuildContinuous\n");
  }

}

void stokesBuildContinuousQuad2D(stokes_t *stokes, dfloat lambda, nonZero_t **A,
				 dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

  mesh_t *mesh = stokes->meshV;
  setupAide options = stokes->options;

  int rank = mesh->rank;

  //use the masked gs handle to define a global ordering

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = stokes->ogs->Ngather;
  dlong Ntotal  = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
  int   *owner     = (int *) calloc(Ngather,sizeof(int));

  // every gathered degree of freedom has its own global id
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<stokes->ogs->Ngather;n++) {
    globalIds[n] = n + globalStarts[rank];
    owner[n] = rank;
  }

  //scatter this numbering to the original nodes
  hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
  for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
  ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, stokes->ogs);
  ogsScatter(globalOwners, owner, ogsInt, ogsAdd, stokes->ogs);

  hlong localMaxGlobalNumber = -100; // ?
  for(dlong n=0;n<Ntotal;++n){
    localMaxGlobalNumber = mymax(localMaxGlobalNumber, globalNumbering[n]);
  }

  // find max global number over all nodes (use this to space stokes entries)
  hlong maxGlobalNumber = -1;
  MPI_Allreduce(&localMaxGlobalNumber, &maxGlobalNumber, 1, MPI_HLONG, MPI_MAX, mesh->comm);
  
  free(globalIds); free(owner);

  // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements*mesh->dim*(1+2);// [ x 0 x; 0 x x; x x 0 ]
  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

  int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
  for (dlong n=0;n<stokes->Nmasked;n++) mask[stokes->maskIds[n]] = 1;

  if(mesh->rank==0) printf("Building full FEM matrix...");fflush(stdout);

  //Build unassembed non-zeros
  // TW: THIS IS ONLY TO BE USED FOR COARSE MESH CONSTRUCTION
  dlong cnt =0;
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int ny=0;ny<mesh->Nq;ny++) {
      for (int nx=0;nx<mesh->Nq;nx++) {

        for (int my=0;my<mesh->Nq;my++) {
          for (int mx=0;mx<mesh->Nq;mx++) {
            
            dfloat Auu = 0, Avv = 0, Aup = 0, Avp = 0, Apu = 0, Apv = 0;

	    for(int ky=0;ky<mesh->Nq;++ky){
	      for(int kx=0;kx<mesh->Nq;++kx){

		int id = kx+ky*mesh->Nq;
		dlong gbase = e*mesh->Np*mesh->Nggeo + id;

		dfloat Grr = mesh->ggeo[gbase + G00ID*mesh->Np];
		dfloat Grs = mesh->ggeo[gbase + G01ID*mesh->Np];
		dfloat Gss = mesh->ggeo[gbase + G11ID*mesh->Np];
		
		dlong vbase = e*mesh->Np*mesh->Nvgeo + id;

		dfloat rx = mesh->vgeo[vbase + RXID*mesh->Np];
		dfloat ry = mesh->vgeo[vbase + RYID*mesh->Np];
		dfloat sx = mesh->vgeo[vbase + SXID*mesh->Np];
		dfloat sy = mesh->vgeo[vbase + SYID*mesh->Np];
		dfloat JW = mesh->vgeo[vbase + JWID*mesh->Np];

		dfloat Bn = (kx==nx)*(ky==ny);
		dfloat Bm = (kx==mx)*(ky==my);
		dfloat Brn = mesh->D[nx+kx*mesh->Nq]*(ky==ny);
		dfloat Bsn = (kx==nx)*mesh->D[ny+ky*mesh->Nq];
		dfloat Brm = mesh->D[mx+kx*mesh->Nq]*(ky==my);
		dfloat Bsm = (kx==mx)*mesh->D[my+ky*mesh->Nq];

		dfloat Bxn = rx*Brn + sx*Bsn;
		dfloat Byn = ry*Brn + sy*Bsn;

		dfloat Bxm = rx*Brm + sx*Bsm;
		dfloat Bym = ry*Brm + sy*Bsm;
		
		Auu += Bxn*JW*Bxm + Byn*JW*Bym; 
		Avv += Bxn*JW*Bxm + Byn*JW*Bym; // clone

		Auu += JW*lambda*Bn*Bm;
		Avv += JW*lambda*Bn*Bm;
		
		Aup += Bxn*JW*Bm; // TW check sign  (dBn/dx, Bm)
		Avp += Byn*JW*Bm; // TW check sign
		
		Apu += Bn*JW*Bxm; // TW check sign
		Apv += Bn*JW*Bym; // TW check sign
		
	      }
	    }

            dfloat nonZeroThreshold = 1e-7;
	    
	    hlong r = globalNumbering[e*mesh->Np + nx+ny*mesh->Nq];
	    hlong c = globalNumbering[e*mesh->Np + mx+my*mesh->Nq];
	    hlong owner = globalOwners[e*mesh->Np + nx+ny*mesh->Nq];
	    
	    // velocity-velocity blocks
	    if (!mask[e*mesh->Np + nx+ny*mesh->Nq]){
	      if (!mask[e*mesh->Np + mx+my*mesh->Nq]){

		if (fabs(Auu)>nonZeroThreshold) {
		  sendNonZeros[cnt].val = Auu;
		  sendNonZeros[cnt].row = r;
		  sendNonZeros[cnt].col = c;
		  sendNonZeros[cnt].ownerRank = owner;	      
		  cnt++;
		}

		if (fabs(Avv)>nonZeroThreshold) {
		  sendNonZeros[cnt].val = Auu; // clone
		  sendNonZeros[cnt].row = r + maxGlobalNumber;
		  sendNonZeros[cnt].col = c + maxGlobalNumber;
		  sendNonZeros[cnt].ownerRank = owner;
		  ++cnt;
		}
	      }
	    }

	    // pressure-velocity blocks
	    if (fabs(Apu)>nonZeroThreshold) {
	      sendNonZeros[cnt].val = Apu;
	      sendNonZeros[cnt].row = r+2*maxGlobalNumber;
	      sendNonZeros[cnt].col = c;
	      sendNonZeros[cnt].ownerRank = owner;	      
	      cnt++;
	    }

	    // pressure-velocity blocks
	    if (fabs(Apv)>nonZeroThreshold) {
	      sendNonZeros[cnt].val = Apv;
	      sendNonZeros[cnt].row = r+2*maxGlobalNumber;
	      sendNonZeros[cnt].col = c+1*maxGlobalNumber;
	      sendNonZeros[cnt].ownerRank = owner;	      
	      cnt++;
	    }

	    // pressure-velocity blocks
	    if (!mask[e*mesh->Np + nx+ny*mesh->Nq]){
	      if (fabs(Aup)>nonZeroThreshold) {
		sendNonZeros[cnt].val = Aup;
		sendNonZeros[cnt].row = r;
		sendNonZeros[cnt].col = c+2*maxGlobalNumber;
		sendNonZeros[cnt].ownerRank = owner;	      
		cnt++;
	      }
	    }
	      
	    // pressure-velocity blocks
	    if (!mask[e*mesh->Np + nx+ny*mesh->Nq]){
	      if (fabs(Avp)>nonZeroThreshold) {
		sendNonZeros[cnt].val = Avp;
		sendNonZeros[cnt].row = r+1*maxGlobalNumber;
		sendNonZeros[cnt].col = c+2*maxGlobalNumber;
		sendNonZeros[cnt].ownerRank = owner;	      
		cnt++;
	      }
	    }
          }
        }
      }
    }
  }
  
  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
  int blength[4] = {1, 1, 1, 1};
  MPI_Aint addr[4], displ[4];
  MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
  MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
  MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
  MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  // count how many non-zeros to send to each process
  for(dlong n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank]++;

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), stokesParallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

  // find send and recv offsets for gather
  *nnz = 0;
  for(int r=0;r<mesh->size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r];
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
		(*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
		mesh->comm);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), stokesParallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<*nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  if (*nnz) cnt++;
  *nnz = cnt;


#if 0
  // Write matlab dat for postprocess 
  char fname[BUFSIZ];
  sprintf(fname, "Ax.dat");
  FILE *fp; 
  fp = fopen(fname, "w");

  for(dlong n=1;n<*nnz;++n){
    fprintf(fp,"%d %d %.8e\n", (*A)[n].row+1, (*A)[n].col+1, (*A)[n].val);
  }

  fclose(fp);
#endif

  if(mesh->rank==0) printf("done.\n");

  MPI_Barrier(mesh->comm);
  MPI_Type_free(&MPI_NONZERO_T);

  free(sendNonZeros);
  free(globalNumbering); free(globalOwners);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
}

