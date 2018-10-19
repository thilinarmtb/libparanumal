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

#include "asbf.h"

int asbfSolve(asbf_t *asbf, setupAide options)
{
  mesh_t *mesh         = asbf->mesh;
  mesh_t *meshSEM      = asbf->meshSEM;
  elliptic_t *elliptic = asbf->pSolver;

  mesh->q = (dfloat*) calloc(mesh->Np*(mesh->Nelements+mesh->totalHaloPairs), sizeof(dfloat));
  for(int m=0;m<asbf->asbfNmodes;++m){
    asbf->o_r.copyFrom(asbf->r3D + asbf->Ntotal*m);

    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      ogsGatherScatter(asbf->o_r, ogsDfloat, ogsAdd, mesh->ogs);
    }

    dfloat lambdam = asbf->lambda + asbf->asbfEigenvalues[m];
    
    // build precon
    if(m>0){
      elliptic->options.setArgs("PRECONDITIONER","JACOBI");
      dfloat *invDiagA;
      ellipticBuildJacobi(elliptic,lambdam,&invDiagA);
      elliptic->precon->o_invDiagA =
	mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
      free(invDiagA);
    }
    
    ellipticSolve(elliptic, lambdam, asbf->pTOL, asbf->o_r, asbf->o_x);
    asbf->o_x.copyTo(asbf->q3D + asbf->Ntotal*m);

    asbf->o_x.copyTo(mesh->q);
    char fileName2D[BUFSIZ];
    sprintf(fileName2D, "bah_%05d.vtu", m);
    meshPlotVTU3D(mesh, fileName2D, 0);
  }

  // reconstruct SEM3D solution ( should do in kernel )
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      // interpolate from asbf to gll nodes
      for(int g=0;g<mesh->Nq;++g){
	dfloat qg = 0;
	for(int i=0;i<asbf->asbfNmodes;++i){
	  qg += asbf->asbfBgll[i + g*asbf->asbfNmodes]*asbf->q3D[(e*mesh->Np+n)+i*asbf->Ntotal];
	}
	// assume Nfields=1
	meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
      }
    }
  }
  
  char fname[] = "foo";
  ellipticPlotVTUHex3D(meshSEM, fname, 0);
}
