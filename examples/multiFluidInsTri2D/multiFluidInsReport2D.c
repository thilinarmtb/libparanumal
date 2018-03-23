#include "multiFluidIns2D.h"

void multiFluidInsReport2D(multiFluidIns_t *solver, int tstep, char *options){

  dfloat t = (tstep+1)*solver->dt;
  
  // copy data back to host
  solver->o_U.copyTo(solver->U);
  solver->o_V.copyTo(solver->V);
  // solver->o_P.copyTo(solver->P);  
  solver->o_Phi.copyTo(solver->Phi);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // // do error stuff on host
  multiFluidInsError2D(solver, t, options);

 
  if(strstr(options, "VTU")){ 
    // compute vorticity
    //insComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins2D/foo_%04d_%04d.vtu",rank, tstep/solver->errorStep);
    sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/solver->errorStep);

    multiFluidInsPlotVTU2D(solver, fname);
  } 
}

