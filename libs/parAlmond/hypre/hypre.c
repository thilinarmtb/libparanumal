#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "omp.h"
#include "crs_hypre.h"

double hypre_param[HYPRE_NPARAM];

struct hypre_crs_data *hypre_setup(int nrows, const long long int rowStart,
                 int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                 const int null_space, const MPI_Comm ce, int Nthreads,
                 const double *param)
{
  struct hypre_crs_data *hypre_data = (struct hypre_crs_data*) malloc(sizeof(struct hypre_crs_data));

  hypre_data->Nthreads = Nthreads;   

  MPI_Comm comm;
  MPI_Comm_dup(ce, &comm);
  int rank;
  MPI_Comm_rank(comm,&rank);

  // Build CRS matrix in hypre format and initialize stuff
  hypre_data->nRows = nrows;
  HYPRE_BigInt ilower = (HYPRE_BigInt) rowStart;
  hypre_data->ilower = ilower;
  HYPRE_BigInt iupper = ilower + (HYPRE_BigInt) nrows - 1; 

  HYPRE_IJMatrixCreate(comm,ilower,iupper,ilower,iupper,&hypre_data->A);
  HYPRE_IJMatrix A_ij = hypre_data->A;
  HYPRE_IJMatrixSetObjectType(A_ij,HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A_ij);

  int i;
  for(i=0; i<nz; i++) 
  {
    HYPRE_BigInt mati = (HYPRE_BigInt)(Ai[i]);
    HYPRE_BigInt matj = (HYPRE_BigInt)(Aj[i]);
    HYPRE_Real matv = (HYPRE_Real) Av[i]; 
    HYPRE_Int ncols = 1; // number of columns per row
    HYPRE_IJMatrixSetValues(A_ij, 1, &ncols, &mati, &matj, &matv);
  }

  HYPRE_IJMatrixAssemble(A_ij);
  //HYPRE_IJMatrixPrint(A_ij, "matrix.dat");

  // Create AMG solver
  HYPRE_BoomerAMGCreate(&hypre_data->solver);
  HYPRE_Solver solver = hypre_data->solver;

  int uparam = (int) param[0];
 
  // Set AMG parameters
  if (uparam) {
      //printf("\n");
      int i;
      for (i = 0; i < HYPRE_NPARAM; i++) {
          hypre_param[i] = param[i+1]; 
           //if (rank == 0)
           //   printf("Custom HYPREsettings[%d]: %.2f\n", i+1, hypre_param[i]);
      }
  } else {
      hypre_param[0] = 10;   /* HMIS                        */
      hypre_param[1] = 6;    /* Extended+i                  */
      hypre_param[2] = 0;    /* not used                    */
      hypre_param[3] = 3;    /* SSOR smoother for crs level */
      hypre_param[4] = 3;
      hypre_param[5] = 0.25;
      hypre_param[6] = 0.1;
      hypre_param[7] = 0.0;
      hypre_param[8] = 0.01;
      hypre_param[9] = 0.05;
  }

  //HYPRE_BoomerAMGSetCycleType(solver,2);
  HYPRE_BoomerAMGSetMaxCoarseSize(solver, 10);
  HYPRE_BoomerAMGSetCoarsenType(solver,hypre_param[0]);
  HYPRE_BoomerAMGSetInterpType(solver,hypre_param[1]);

  HYPRE_BoomerAMGSetKeepTranspose(solver, 1);

  HYPRE_BoomerAMGSetCycleRelaxType(solver, 16, 1);
  HYPRE_BoomerAMGSetCycleRelaxType(solver, 16, 2);
  HYPRE_BoomerAMGSetCycleRelaxType(solver,  9, 3);
  HYPRE_BoomerAMGSetChebyFraction(solver, 1./20); 

  if (null_space) {
    HYPRE_BoomerAMGSetCycleRelaxType(solver, hypre_param[3], 3);
    HYPRE_BoomerAMGSetCycleNumSweeps(solver, hypre_param[4], 3);
  }

  HYPRE_BoomerAMGSetStrongThreshold(solver,hypre_param[5]);

  HYPRE_BoomerAMGSetNonGalerkinTol(solver,hypre_param[6]);
  HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,hypre_param[7], 0);
  HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,hypre_param[8], 1);
  HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,hypre_param[9], 2);
 
  //HYPRE_BoomerAMGSetAggNumLevels(solver, 1); 

  HYPRE_BoomerAMGSetMaxIter(solver,hypre_param[2]); // number of V-cycles
  HYPRE_BoomerAMGSetTol(solver,0);

  HYPRE_BoomerAMGSetPrintLevel(solver,1);

  // Create and initialize rhs and solution vectors
  HYPRE_IJVectorCreate(comm,ilower,iupper,&hypre_data->b);
  HYPRE_IJVector b = hypre_data->b;
  HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);
  HYPRE_IJVectorAssemble(b);

  HYPRE_IJVectorCreate(comm,ilower,iupper,&hypre_data->x);
  HYPRE_IJVector x = hypre_data->x;
  HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);
  HYPRE_IJVectorAssemble(x);

  // Perform AMG setup
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;
  HYPRE_IJVectorGetObject(b,(void**) &par_b);
  HYPRE_IJVectorGetObject(x,(void**) &par_x);
  HYPRE_ParCSRMatrix par_A;
  HYPRE_IJMatrixGetObject(hypre_data->A,(void**) &par_A);

  int _Nthreads = omp_get_thread_num();
  omp_set_num_threads(hypre_data->Nthreads);
  HYPRE_BoomerAMGSetup(solver,par_A,par_b,par_x);
  omp_set_num_threads(_Nthreads);

  hypre_data->ii = (HYPRE_BigInt*) malloc(hypre_data->nRows*sizeof(HYPRE_BigInt));
  hypre_data->bb = (HYPRE_Real*) malloc(hypre_data->nRows*sizeof(HYPRE_Real));
  hypre_data->xx = (HYPRE_Real*) malloc(hypre_data->nRows*sizeof(HYPRE_Real));
  for(i=0;i<hypre_data->nRows;++i) 
    hypre_data->ii[i] = ilower + (HYPRE_BigInt)i;

  return hypre_data;
}

void hypre_solve(double *x, struct hypre_crs_data *data, double *b)
{
  int i; 
  HYPRE_BigInt ilower = data->ilower;
  HYPRE_IJVector ij_x = data->x;
  HYPRE_IJVector ij_b = data->b;
  HYPRE_IJMatrix ij_A = data->A;
  HYPRE_Solver solver = data->solver;

  HYPRE_ParVector par_x;
  HYPRE_ParVector par_b;
  HYPRE_ParCSRMatrix par_A;

  for(i=0;i<data->nRows;++i) 
    data->bb[i] = b[i]; 
  HYPRE_IJVectorSetValues(ij_b,data->nRows,data->ii,data->bb);

  HYPRE_IJVectorAssemble(ij_b);
  HYPRE_IJVectorGetObject(ij_b,(void**) &par_b);

  HYPRE_IJVectorAssemble(ij_x);
  HYPRE_IJVectorGetObject(ij_x,(void **) &par_x);

  HYPRE_IJMatrixGetObject(ij_A,(void**) &par_A);

  int _Nthreads = omp_get_thread_num();
  omp_set_num_threads(data->Nthreads);
  HYPRE_BoomerAMGSolve(solver,par_A,par_b,par_x);
  omp_set_num_threads(_Nthreads);

  HYPRE_IJVectorGetValues(ij_x,data->nRows,data->ii,data->xx);
  for(i=0;i<data->nRows;++i) 
    x[i] = data->xx[i]; 
}

void hypre_free(struct hypre_crs_data *data)
{
  HYPRE_BoomerAMGDestroy(data->solver);
  HYPRE_IJMatrixDestroy(data->A);
  HYPRE_IJVectorDestroy(data->x);
  HYPRE_IJVectorDestroy(data->b);
  free(data);
}