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

// compare on global indices
static int parallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

void ellipticBuildContinuousGalerkinHex3D (elliptic_t *elliptic, elliptic_t *ellipticFine,
  dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);

void ellipticBuildContinuousGalerkin(elliptic_t *elliptic, elliptic_t *ellipticFine,
  dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

  mesh_t *mesh=elliptic->mesh;
  switch(elliptic->elementType){
    case TRIANGLES:
    case TETRAHEDRA:
    case QUADRILATERALS:
      printf("ellipticBuildContinuousGalerkin is called with wrong element type.\n");
      exit(1);
      break;
    case HEXAHEDRA:
      ellipticBuildContinuousGalerkinHex3D(elliptic,ellipticFine,lambda,A,nnz,ogs,globalStarts);
      break;

      break;
    default:
      break;
  }
}

#define p_Nggeo 7

template < const int rowsA, const int rowsB, const int colsC >
  static void mxm(const dfloat * __restrict__ A,
		    const dfloat * __restrict__ B,
		    const dfloat BETA, 
		    dfloat * __restrict__ C)
{
  if(BETA)
    for(int j=0;j<colsC;++j){
      for(int i=0;i<rowsA;++i){
	dfloat res = BETA*C[i+j*rowsA];
	for(int k=0;k<rowsB;++k){
	  res += A[i+k*rowsA]*B[k+j*rowsB];
	}
	C[i+j*rowsA] = res;
      }
    }
  else
    for(int j=0;j<colsC;++j){
      for(int i=0;i<rowsA;++i){
	dfloat res = 0;
	for(int k=0;k<rowsB;++k){
	  res += A[i+k*rowsA]*B[k+j*rowsB];
	}
	C[i+j*rowsA] = res;
      }
    }
}

template < const int p_Nq >
void ellipticHostElementAxHexKernel3D(const dfloat * __restrict__ ggeo,
					const dfloat * __restrict__ D,
					const dfloat * __restrict__ S,
					const dfloat * __restrict__ MM,
					const dfloat lambda,
					const dfloat * __restrict__ q,
					dfloat * __restrict__ qr,
					dfloat * __restrict__ qs,
					dfloat * __restrict__ qt,
					dfloat * __restrict__ Aq,
					dfloat * __restrict__ wk)
{
  const int p_Np=p_Nq*p_Nq*p_Nq;
  const int p_N =p_Nq-1;

  D    = (dfloat*)__builtin_assume_aligned(D, USE_OCCA_MEM_BYTE_ALIGN) ;
  S    = (dfloat*)__builtin_assume_aligned(S, USE_OCCA_MEM_BYTE_ALIGN) ;
  MM   = (dfloat*)__builtin_assume_aligned(MM, USE_OCCA_MEM_BYTE_ALIGN) ;
  q    = (dfloat*)__builtin_assume_aligned(q, USE_OCCA_MEM_BYTE_ALIGN) ;
  qr    = (dfloat*)__builtin_assume_aligned(qr, USE_OCCA_MEM_BYTE_ALIGN) ;
  qs    = (dfloat*)__builtin_assume_aligned(qs, USE_OCCA_MEM_BYTE_ALIGN) ;
  qt    = (dfloat*)__builtin_assume_aligned(qt, USE_OCCA_MEM_BYTE_ALIGN) ;
  Aq   = (dfloat*)__builtin_assume_aligned(Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
  ggeo = (dfloat*)__builtin_assume_aligned(ggeo, USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat zero = 0, one = 1.0;
  const int N = p_Nq-1;

  // grad
  mxm<p_Nq,p_Nq,p_Nq*p_Nq>(S, q, zero, qr); // D(:,:)*q(:,:+::)
  for(int k=0;k<p_Nq;++k){
    mxm<p_Nq,p_Nq,p_Nq>(q+k*p_Nq*p_Nq, D, zero, qs+k*p_Nq*p_Nq);
  }
  mxm<p_Nq*p_Nq,p_Nq,p_Nq>(q, D, zero, qt);
  
  for(int n=0;n<p_Np;++n){
    const dfloat G00 = ggeo[n+G00ID*p_Np], G01 = ggeo[n+G01ID*p_Np], G11 = ggeo[n+G11ID*p_Np];
    const dfloat GWJ = ggeo[n+GWJID*p_Np];
    const dfloat G12 = ggeo[n+G12ID*p_Np], G02 = ggeo[n+G02ID*p_Np], G22 = ggeo[n+G22ID*p_Np];

    dfloat qrn = G00*qr[n] + G01*qs[n] + G02*qt[n];
    dfloat qsn = G01*qr[n] + G11*qs[n] + G12*qt[n];
    dfloat qtn = G02*qr[n] + G12*qs[n] + G22*qt[n];

    qr[n] = qrn;
    qs[n] = qsn;
    qt[n] = qtn;
  }

  // gradT
  mxm<p_Nq,p_Nq,p_Nq*p_Nq>(D, qr, zero, Aq); // D(:,:)*q(:,:+::)
  for(int k=0;k<p_Nq;++k){
    mxm<p_Nq,p_Nq,p_Nq>(qs+k*p_Nq*p_Nq, S, one, Aq+k*p_Nq*p_Nq);
  }
  mxm<p_Nq*p_Nq,p_Nq,p_Nq>(qt, S, one, Aq);
}


template < const int p_Nq >
void ellipticHostAxHexKernel3D (const hlong Nelements,
				  const dfloat * __restrict__ ggeo ,
				  const dfloat * __restrict__ D ,
				  const dfloat * __restrict__ S ,
				  const dfloat * __restrict__ MM ,
				  const dfloat lambda,
				  const dfloat * __restrict__ q ,
				  dfloat * __restrict__ Aq )
{
  const int p_Np=p_Nq*p_Nq*p_Nq;
  const int c_Np=p_Np;
  const int p_N =p_Nq-1;

  D    = (dfloat*)__builtin_assume_aligned(D, USE_OCCA_MEM_BYTE_ALIGN) ;
  S    = (dfloat*)__builtin_assume_aligned(S, USE_OCCA_MEM_BYTE_ALIGN) ;
  MM   = (dfloat*)__builtin_assume_aligned(MM, USE_OCCA_MEM_BYTE_ALIGN) ;
  q    = (dfloat*)__builtin_assume_aligned(q, USE_OCCA_MEM_BYTE_ALIGN) ;
  Aq   = (dfloat*)__builtin_assume_aligned(Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
  ggeo = (dfloat*)__builtin_assume_aligned(ggeo, USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat s_tmp[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  
  dfloat s_qr[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_qs[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_qt[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_wk[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));

  for(dlong e=0; e<Nelements; ++e){
    const dlong element = e;

    // MIXED STEFAN + TW VERSION
    ellipticHostElementAxHexKernel3D<p_Nq>(ggeo+element*p_Np*p_Nggeo,
					     D, S, MM, lambda, q + element*p_Np,
					     s_qr, s_qs, s_qt, Aq+element*p_Np, s_wk);
  }
}

void ellipticHostAxHexKernel3D(const int Nq,
				 const hlong Nelements,
				 const dfloat *ggeo,
				 const dfloat *D,
				 const dfloat *S,
				 const dfloat *MM,
				 const dfloat lambda,
				 const dfloat *q,
				       dfloat *Aq
				 )
{
  switch(Nq){
  case  2:
    ellipticHostAxHexKernel3D< 2> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  3:
    ellipticHostAxHexKernel3D< 3> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  4:
    ellipticHostAxHexKernel3D< 4> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  5:
    ellipticHostAxHexKernel3D< 5> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  6:
    ellipticHostAxHexKernel3D< 6> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  7:
    ellipticHostAxHexKernel3D< 7> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  8:
    ellipticHostAxHexKernel3D< 8> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case  9:
    ellipticHostAxHexKernel3D< 9> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case 10:
    ellipticHostAxHexKernel3D<10> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case 11:
    ellipticHostAxHexKernel3D<11> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case 12:
    ellipticHostAxHexKernel3D<12> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  case 13:
    ellipticHostAxHexKernel3D<13> (Nelements,ggeo,D,S,MM,lambda,q,Aq); break;
  }
}

void ellipticGenerateCoarseBasisHex3D(dfloat *b,int j_,elliptic_t *elliptic){
  // b = lx1**dim
  mesh_t* mesh=elliptic->mesh;
  dfloat *z=mesh->gllz;

  int jj=j_+1;
  dfloat *zr,*zs,*zt,*z0,*z1;
  zr=(dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  zs=(dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  zt=(dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  z0=(dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  z1=(dfloat*)calloc(mesh->Nq,sizeof(dfloat));

  for(int i=0;i<mesh->Nq;i++){
    z0[i]=0.5*(1-z[i]);
    z1[i]=0.5*(1+z[i]);
  }

  memcpy(zr,z0,mesh->Nq*sizeof(dfloat));
  memcpy(zs,z0,mesh->Nq*sizeof(dfloat));
  memcpy(zt,z0,mesh->Nq*sizeof(dfloat));

  if(jj%2 == 0)                        memcpy(zr,z1,mesh->Nq*sizeof(dfloat));
  if(jj==3 || jj==4 || jj==7 || jj==8) memcpy(zs,z1,mesh->Nq*sizeof(dfloat));
  if(jj>4)                             memcpy(zt,z1,mesh->Nq*sizeof(dfloat));

  int dim=mesh->dim;
  for(int k=0; k<mesh->Nq; k++){
  for(int j=0; j<mesh->Nq; j++){
  for(int i=0; i<mesh->Nq; i++){
    int n=i+mesh->Nq*j+mesh->Nq*mesh->Nq*k+j_*mesh->Np;
    b[n]=zr[i]*zs[j]*zt[k];
  }
  }
  }
  free(zr); free(zs); free(zt); free(z0); free(z1);
}

void ellipticBuildContinuousGalerkinHex3D(elliptic_t *elliptic,elliptic_t *ellipticFine,
  dfloat lambda,nonZero_t **A,dlong *nnz,ogs_t **ogs,hlong *globalStarts)
{

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rank = mesh->rank;

  //use the masked gs handle to define a global ordering

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = elliptic->ogs->Ngather;
  dlong Ntotal  = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
  int   *owner     = (int *) calloc(Ngather,sizeof(int));

  // every gathered degree of freedom has its own global id
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<elliptic->ogs->Ngather;n++) {
    globalIds[n] = n + globalStarts[rank];
    owner[n] = rank;
  }

  //scatter this numbering to the original nodes
  hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
  for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
  ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, elliptic->ogs);
  ogsScatter(globalOwners, owner, ogsInt, ogsAdd, elliptic->ogs);

  free(globalIds); free(owner);


    // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;
  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

  int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
  for (dlong n=0;n<elliptic->Nmasked;n++) mask[elliptic->maskIds[n]] = 1;

  if(mesh->rank==0) printf("Building full FEM matrix...\n"
    "Using coarse system generated by fine level...\n");
  fflush(stdout);

  mesh_t *meshf=ellipticFine->mesh;

  dfloat *b,*q,*Aq;
  b =(dfloat *)calloc(mesh->Np*meshf->Np       ,sizeof(dfloat));
  q =(dfloat *)calloc(meshf->Np*mesh->Nelements,sizeof(dfloat));
  Aq=(dfloat *)calloc(meshf->Np*mesh->Nelements,sizeof(dfloat));

  for(int jj=0;jj<mesh->Np;jj++)
    ellipticGenerateCoarseBasisHex3D(b,jj,ellipticFine);

  dlong cnt =0;
  for (int nz=0;nz<mesh->Nq;nz++) {
  for (int ny=0;ny<mesh->Nq;ny++) {
  for (int nx=0;nx<mesh->Nq;nx++) {
    int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
    // ok
    for (dlong e=0;e<mesh->Nelements;e++) {
      memcpy(&q[e*meshf->Np],&b[idn*meshf->Np],meshf->Np*sizeof(dfloat));
    }

    // Aq
    ellipticHostAxHexKernel3D(meshf->Nq,meshf->Nelements,meshf->ggeo,meshf->D,
	    meshf->DT,meshf->MM,lambda,q,Aq);

    for(dlong e=0;e<mesh->Nelements;e++){
      for (int mz=0;mz<mesh->Nq;mz++) {
      for (int my=0;my<mesh->Nq;my++) {
      for (int mx=0;mx<mesh->Nq;mx++) {
         int idm = mx+my*mesh->Nq+mz*mesh->Nq*mesh->Nq;
         if (mask[e*mesh->Np + idm])
           continue; //skip masked nodes
         if (mask[e*mesh->Np + idn])
           continue; //skip masked nodes

         dfloat val=0.f;
         for(int mmm=0;mmm<meshf->Np;mmm++){
           val+=Aq[e*meshf->Np+mmm]*b[mmm+idm*meshf->Np];
         }

         // pack non-zero
         dfloat nonZeroThreshold = 1e-7;
         if (fabs(val) >= nonZeroThreshold) {
           sendNonZeros[cnt].val = val;
           sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + idm];
           sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + idn];
           sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Np + idm];
           cnt++;
         }
      }
      }
      }
    }
  }
  }
  }

  free(b); free(q); free(Aq);

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
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

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
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

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
