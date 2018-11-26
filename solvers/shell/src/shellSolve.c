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

int shellSolve(shell_t *shell, setupAide options)
{
  mesh_t *mesh         = shell->mesh;
  mesh_t *meshSEM      = shell->meshSEM;
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

    // build precon
    if(m>5){
      elliptic->options.setArgs("PRECONDITIONER","JACOBI");
      dfloat *invDiagA;
      ellipticBuildJacobi(elliptic,lambdam,&invDiagA);
      elliptic->precon->o_invDiagA =
        mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
      free(invDiagA);
    }

    ellipticSolve(elliptic, lambdam, shell->pTOL, shell->o_r, shell->o_x);
    shell->o_x.copyTo(shell->q3D + shell->Ntotal*m);

    shell->o_x.copyTo(mesh->q);

    // Plot solutions for each mode.
    //char fileName2D[BUFSIZ];
    //sprintf(fileName2D, "bah_%05d.vtu", m);
    //meshPlotVTU3D(mesh, fileName2D, 0);
  }

#if 0
  // reconstruct SEM3D solution ( should do in kernel )
  if(shell->elementType==QUADRILATERALS){
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        // interpolate from shell to gll nodes
        for(int g=0;g<mesh->Nq;++g){
          dfloat qg = 0;
          for(int i=0;i<shell->Nmodes;++i){
            qg += shell->Bgll[i + g*shell->Nmodes]*shell->q3D[(e*mesh->Np+n)+i*shell->Ntotal];
          }
          // assume Nfields=1
          meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
        }
      }
    }
  }
#endif  
}
