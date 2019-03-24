
#include "stokes.h"

/*
To test:

make sparseMxM

 ./sparseMxM 8000 8000 200000 8000 8000 200000

in matlab:

A = load('A00000.dat');
B = load('B00000.dat');
C = load('C00000.dat');

A = spconvert(A);
B = spconvert(B);
C = spconvert(C);

max(max(abs(A*B-C)))


 */

typedef struct {
  
  int matrixId; // 0 = A, 1 = B, 2 = C
  hlong  row;
  hlong  col;
  dfloat val;
  hlong sortId;
  int   source;

}entry_t;

int compareSortId(const void *a, const void *b){

  entry_t *eA = (entry_t*) a;
  entry_t *eB = (entry_t*) b;

  if(eA->sortId < eB->sortId) return -1;
  if(eA->sortId > eB->sortId) return +1;
  return 0;
}

int compareIndices(const void *a, const void *b){

  entry_t *eA = (entry_t*) a;
  entry_t *eB = (entry_t*) b;

  if(eA->matrixId < eB->matrixId) return -1;
  if(eA->matrixId > eB->matrixId) return +1;
  
  if(eA->row < eB->row) return -1;
  if(eA->row > eB->row) return +1;

  if(eA->col < eB->col) return -1;
  if(eA->col > eB->col) return +1;

  return 0;
}

// ------------------------------------------------------------------------------------
// accumulate duplicates
hlong sparseAccumulateDuplicates(hlong Nentries, entry_t *entries){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  qsort(entries, Nentries, sizeof(entry_t), compareIndices); // sort by col, then row

  hlong cnt = 0;
  hlong n = 0, m;
  while(n<Nentries){
    dfloat val = 0;
    for(m=n;m<Nentries;++m){
      
      if(entries[n].matrixId == entries[m].matrixId &&
	 entries[n].row == entries[m].row &&
	 entries[n].col == entries[m].col){
	val += entries[m].val;
      }
      else{
	break;
      }
    }
    entries[cnt] = entries[n];
    entries[cnt].val = val;
    n = m;
    ++cnt;
  }
  
  return cnt;
}


hlong sparseRedistribute(hlong NsendEntries,
			 entry_t *sendEntries,
			 entry_t **recvEntries){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  hlong n;

  hlong *sendCounts = (hlong*) calloc(size, sizeof(hlong));
  hlong *recvCounts = (hlong*) calloc(size, sizeof(hlong));

  hlong *sendOffsets = (hlong*) calloc(size+1, sizeof(hlong));
  hlong *recvOffsets = (hlong*) calloc(size+1, sizeof(hlong));

  qsort(sendEntries, NsendEntries, sizeof(entry_t), compareSortId); // sort by sortId
  
  for(n=0;n<NsendEntries;++n){
    ++sendCounts[sendEntries[n].sortId];
  }

  MPI_Alltoall(sendCounts, 1, MPI_HLONG,
	       recvCounts, 1, MPI_HLONG, MPI_COMM_WORLD);

  for(n=0;n<size;++n){

    sendCounts[n] *= sizeof(entry_t);
    recvCounts[n] *= sizeof(entry_t);

    sendOffsets[n+1] = sendOffsets[n] + sendCounts[n];
    recvOffsets[n+1] = recvOffsets[n] + recvCounts[n];

  }

  hlong cnt = recvOffsets[size]/sizeof(entry_t);
  *recvEntries = (entry_t*) calloc(cnt, sizeof(entry_t));

  MPI_Alltoallv( sendEntries, sendCounts, sendOffsets, MPI_CHAR,
		*recvEntries, recvCounts, recvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  free(sendCounts);
  free(recvCounts);
  free(sendOffsets);
  free(recvOffsets);

  return cnt;

}



int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  // demo code for distributed mxm
  hlong NrowsA = atoi(argv[1]);
  hlong NcolsA = atoi(argv[2]);
  hlong NrandA = atoi(argv[3]); // nnz per rank
  
  hlong NrowsB = atoi(argv[4]);
  hlong NcolsB = atoi(argv[5]);
  hlong NrandB = atoi(argv[6]); // nnz per rank
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // seed change on rank
  srand48(rank*12345);
  
  // make no assumptions about distribution
  dfloat *Avals = (dfloat*) calloc(NrandA, sizeof(dfloat));
  hlong  *Arows = (hlong*)  calloc(NrandA, sizeof(hlong));
  hlong  *Acols = (hlong*)  calloc(NrandA, sizeof(hlong)); // ok - triplet, not csr

  dfloat *Bvals = (dfloat*) calloc(NrandB, sizeof(dfloat));
  hlong  *Brows = (hlong*)  calloc(NrandB, sizeof(hlong));
  hlong  *Bcols = (hlong*)  calloc(NrandB, sizeof(hlong)); // ok - triplet, not csr

  // populate A (local entries)
  for(hlong n=0;n<NrandA;++n){
    Arows[n] = mymin((hlong)(drand48()*NrowsA), NrowsA-1);
    Acols[n] = mymin((hlong)(drand48()*NcolsA), NcolsA-1);
    Avals[n] = drand48();
  }

  // populate B (local entries)
  for(hlong n=0;n<NrandB;++n){
    Brows[n] = mymin((hlong)(drand48()*NrowsB), NrowsB-1);
    Bcols[n] = mymin((hlong)(drand48()*NcolsB), NcolsB-1);
    Bvals[n] = drand48();
  }
  
  // GOAL: compute C = A*B
  
  // 1. merge A and B into combined array
  entry_t *sendEntries = (entry_t*) calloc(NrandA+NrandB, sizeof(entry_t));
  hlong cnt = 0;
  for(hlong n=0;n<NrandA;++n){
    sendEntries[cnt].matrixId = 0;   // A label
    sendEntries[cnt].row = Arows[n];
    sendEntries[cnt].col = Acols[n];
    sendEntries[cnt].val = Avals[n];
    sendEntries[cnt].source = rank;
    sendEntries[cnt].sortId = sendEntries[cnt].col%size; // make sure entries from same column of A sent to same rank
    ++cnt;
  }

  for(hlong n=0;n<NrandB;++n){
    sendEntries[cnt].matrixId = 1;   // B label
    sendEntries[cnt].row = Brows[n]; 
    sendEntries[cnt].col = Bcols[n];
    sendEntries[cnt].val = Bvals[n];
    sendEntries[cnt].source = rank;
    sendEntries[cnt].sortId = sendEntries[cnt].row%size; // make sure entries from same row of B sent to same rank
    ++cnt;
  }

  // ------------------------------------------------------------------------------------
  entry_t *recvEntries;

  // redistribute by sortId (A col index, or B row index)
  cnt = sparseRedistribute(cnt, sendEntries, &recvEntries);

  // remove duplicates (complete cols of A and rows of B now on same ranks)
  cnt = sparseAccumulateDuplicates(cnt, recvEntries);

  // ------------------------------------------------------------------------------------
  // resort to make sure all entries with the same column (A) and row (B) are adjacent
  for(hlong n=0;n<cnt;++n){
    if(recvEntries[n].matrixId==0)  // A
      recvEntries[n].sortId = recvEntries[n].col;
    if(recvEntries[n].matrixId==1)  // B
      recvEntries[n].sortId = recvEntries[n].row;
  }

  qsort(recvEntries, cnt, sizeof(entry_t), compareSortId); // sort by sortId

  // ------------------------------------------------------------------------------------
  // count number of values produced in matrix-matrix multiplication
  hlong newCount = 0, n=0, m;
  while(n<cnt){
    
    for(m=n+1;m<cnt;++m){
      if(m<cnt){
	if(recvEntries[m].sortId!=recvEntries[n].sortId){
	  break;
	}
      }
      else
	break;
    }
    
    // loop through all entries with the same col 
    hlong start = n;
    hlong end = m;
    
    // find pair of A and B entries
    for(hlong r=start;r<end;++r){
      if(recvEntries[r].matrixId == 0) {
	for(hlong c=start;c<end;++c){
	  if(recvEntries[c].matrixId == 1) {
	    ++newCount;
	  }
	}
      }
    }
    
    n = m;
  }
  
  // ------------------------------------------------------------------------------------
  // Perform matrix-matrix product without accumulating
  // TW: leak
  sendEntries = (entry_t*) calloc(newCount,sizeof(entry_t));

  // reset counters
  newCount = 0;
  n = 0;
  while(n<cnt){
    
    for(m=n+1;m<cnt;++m){
      if(m<cnt){
	if(recvEntries[m].sortId!=recvEntries[n].sortId){
	  break;
	}
      }
      else
	break;
    }
    
    // loop through all entries with the same col 
    hlong start = n;
    hlong end = m;

    // find pair of A and B entries
    for(hlong r=start;r<end;++r){
      if(recvEntries[r].matrixId == 0) {
	for(hlong c=start;c<end;++c){
	  if(recvEntries[c].matrixId == 1) {
	    
	    sendEntries[newCount].matrixId = 2; // C
	    sendEntries[newCount].row = recvEntries[r].row;
	    sendEntries[newCount].col = recvEntries[c].col; 
	    sendEntries[newCount].val = recvEntries[r].val*recvEntries[c].val;
	    sendEntries[newCount].source = recvEntries[r].source;
	    ++newCount;
	  }
	}
      }
    }
    
    n = m;
  }


  // ------------------------------------------------------------------------------------
  // accumulate duplicates
  newCount = sparseAccumulateDuplicates(newCount, sendEntries);

  // ------------------------------------------------------------------------------------
  // send non-zero entries to intermediate for accumulation (distribute by sortId)

  for(n=0;n<newCount;++n){
    sendEntries[n].sortId = sendEntries[n].row%size; 
  }
  
  newCount = sparseRedistribute(newCount, sendEntries, &recvEntries);

  // ------------------------------------------------------------------------------------
  // accumulate duplicates
  newCount = sparseAccumulateDuplicates(newCount, recvEntries);


  char nameA[BUFSIZ];
  char nameB[BUFSIZ];
  char nameC[BUFSIZ];

  sprintf(nameA, "A%05d.dat", rank);
  sprintf(nameB, "B%05d.dat", rank);
  sprintf(nameC, "C%05d.dat", rank);
    
  FILE *fpA = fopen(nameA, "w");
  FILE *fpB = fopen(nameB, "w");
  FILE *fpC = fopen(nameC, "w");

  for(n=0;n<NrandA;++n){
    fprintf(fpA, hlongFormat " " hlongFormat " %17.15lg \n",
	    Arows[n]+1,
	    Acols[n]+1,
	    Avals[n]);
  }

  for(n=0;n<NrandB;++n){
    fprintf(fpB, hlongFormat " " hlongFormat " %17.15lg \n",
	    Brows[n]+1,
	    Bcols[n]+1,
	    Bvals[n]);
  }

  for(n=0;n<newCount;++n){
    fprintf(fpC, hlongFormat " " hlongFormat " %17.15lg \n",
	    recvEntries[n].row+1,
	    recvEntries[n].col+1,
	    recvEntries[n].val);
  }
  
  fclose(fpA);
  fclose(fpB);
  fclose(fpC); 
  
  MPI_Finalize();
  exit(0);
  return 0;
  
}
