#ifndef CRS_HYPRE_H
#define CRS_HYPRE_H

#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_ls.h"
#include "HYPRE.h"

#define HYPRE_NPARAM 10

struct hypre_crs_data {
  HYPRE_Solver solver;
  HYPRE_IJMatrix A;
  HYPRE_IJVector b;
  HYPRE_IJVector x;
  HYPRE_BigInt ilower;
  HYPRE_BigInt *ii;
  HYPRE_Real *bb;
  HYPRE_Real *xx;
  int nRows;
  int Nthreads; 
};

#ifdef __cplusplus
extern "C" {
#endif

struct hypre_crs_data *hypre_setup(int nrows, const long long int rowStart,
                 int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                 const int null_space, const MPI_Comm ce, int Nthreads,
                 const double *param);

void hypre_solve(double *x, struct hypre_crs_data *data, double *b);

void hypre_free(struct hypre_crs_data *data);

#ifdef __cplusplus
}
#endif

#endif
