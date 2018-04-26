#include "cnsQuad2D.h"

void cnsReportQuad2D(cns_t *cns, int tstep, char* options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh2D *mesh = cns->mesh;

  dfloat t = (tstep)*mesh->dt;
  
  cns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_D,
                       cns->o_q,
                       cns->o_Vort);

  // copy data back to host
  cns->o_q.copyTo(mesh->q);
  cns->o_Vort.copyTo(cns->Vort);

  // do error stuff on host
  cnsError2D(mesh, mesh->dt*(tstep+1));


  // output field files
  char fname[BUFSIZ];
  // sprintf(fname, "/u0/outputs/cns2D/foo_%04d_%04d.vtu",rank, tstep/cns->errorStep);
  sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/mesh->errorStep);
  //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank, tstep/cns->errorStep);
  cnsPlotVTUQuad2D(cns, fname);

}