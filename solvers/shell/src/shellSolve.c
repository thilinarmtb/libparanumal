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

static int shellSolveASBF(shell_t *shell, setupAide options);
static int shellSolveSEM(shell_t *shell, setupAide options);

int shellSolve(shell_t *shell, setupAide options)
{
  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    return shellSolveASBF(shell, options);
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    return shellSolveSEM(shell, options);
  } else {
    printf("ERROR:  Invalid value \"%s\" for SHELL SOLVER.\n",
           options.getArgs("SHELL SOLVER").c_str());
    exit(-1);
  }
}

static int shellSolveASBF(shell_t *shell, setupAide options)
{
  mesh_t *mesh         = shell->mesh;
  elliptic_t *elliptic = shell->elliptic;

  mesh->q = (dfloat*) calloc(mesh->Np*(mesh->Nelements+mesh->totalHaloPairs), sizeof(dfloat));
  for(int m=0;m<shell->Nmodes;++m){

    // integrate agains surface basis (assume scaled by J)
    if (options.compareArgs("BASIS","NODAL")){
      meshApplyElementMatrix(mesh,mesh->MM,shell->r3D+shell->Ntotal*m, shell->r3D+shell->Ntotal*m);
    }

    shell->o_r.copyFrom(shell->r3D + shell->Ntotal*m);

    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      // gather scatter
      ogsGatherScatter(shell->o_r, ogsDfloat, ogsAdd, mesh->ogs);
    }

    //dfloat lambdam = shell->lambda + shell->eigenvalues[m];
    dfloat lambdam = shell->eigenvalues[m];
    printf("LAMBDAM[%02d] = %.3f\n", m, lambdam);

    // Switch to JACOBI for high-order modes.
    if (m >= SHELL_ASBF_JACOBI_CROSSOVER) {
      elliptic->options.setArgs("PRECONDITIONER", "JACOBI");
      elliptic->precon = shell->preconJacobi + (m - SHELL_ASBF_JACOBI_CROSSOVER);
    }

    ellipticSolve(elliptic, lambdam, shell->TOL, shell->o_r, shell->o_x);
    shell->o_x.copyTo(shell->q3D + shell->Ntotal*m);

    shell->o_x.copyTo(mesh->q);

    // Plot solutions for each mode.
    //char fileName2D[BUFSIZ];
    //sprintf(fileName2D, "bah_%05d.vtu", m);
    //meshPlotVTU3D(mesh, fileName2D, 0);
  }
}

static int shellSolveSEM(shell_t *shell, setupAide options)
{
  int Niters;

  elliptic_t *elliptic = shell->elliptic;
  mesh_t *mesh         = elliptic->mesh;

  if (options.compareArgs("BASIS", "NODAL"))
    meshApplyElementMatrix(mesh, mesh->MM, elliptic->r, elliptic->r);

  elliptic->o_r.copyFrom(elliptic->r);

  if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    ogsGatherScatter(elliptic->o_r, ogsDfloat, ogsAdd, mesh->ogs);
    if (elliptic->Nmasked)
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_r);
  }
  
  Niters = ellipticSolve(elliptic, shell->lambda, shell->TOL, elliptic->o_r, elliptic->o_x);
  elliptic->o_x.copyTo(mesh->q);

  return Niters;
}
