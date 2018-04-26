#include "cnsQuad2D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // SET OPTIONS
  // out  = REPORT, REPORT+VTU
  // adv  = CUBATURE, COLLOCATION
  char *options = strdup("out=VTU, adv=CUBATURE"); 

  // set up mesh stuff
  mesh2D *mesh = meshSetupQuad2D(argv[1], N);

  char *boundaryHeaderFileName = strdup(DHOLMES "/examples/cnsQuad2D/cnsUniform2D.h"); // default

  // set up cns stuff
  cns_t *cns = cnsSetupQuad2D(mesh, options, boundaryHeaderFileName);

  // run
  cnsRunQuad2D(cns, options);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}