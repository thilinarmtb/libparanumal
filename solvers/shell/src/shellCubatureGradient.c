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

// TODO:  Doesn't this belong in, e.g., meshGeometricFactorsQuad2D()?
void interpolateQuad2D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M)
{
  dfloat *Ix1 = (dfloat*) calloc(N*M, sizeof(dfloat));

  for(int j=0;j<N;++j){
    for(int i=0;i<M;++i){
      dfloat tmp = 0;
      for(int n=0;n<N;++n){
        tmp += I[i*N + n]*x[j*N+n];
      }
      Ix1[j*M+i] = tmp;
    }
  }

  for(int j=0;j<M;++j){
    for(int i=0;i<M;++i){
      dfloat tmp = 0;
      for(int n=0;n<N;++n){
        tmp += I[j*N + n]*Ix1[n*M+i];
      }
      Ix[j*M+i] = tmp;
    }
  }

  free(Ix1);
}

void shellCubatureGradient(shell_t *shell, dfloat *q3D,
    dfloat *cubq, dfloat *cubdqdx, dfloat *cubdqdy, dfloat *cubdqdz){

  mesh_t *mesh = shell->mesh;

  dfloat *qe  = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *cubqe  = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqre = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqse = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqte = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  dfloat* cubx = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cuby = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubz = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  for(int e = 0; e < mesh->Nelements; e++){

    // interpolate and differentiate
    for(int k=0;k<shell->Nquad;++k){

      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){

          int m = j*mesh->Nq + i;
          dfloat qk = 0, dqk = 0;
          for(int p=0;p<shell->Nmodes;++p){
            dfloat qpji = shell->q3D[(e*mesh->Np+m)+p*shell->Ntotal];
            qk  +=  shell->Bquad[p+k*shell->Nmodes]*qpji;
            dqk += shell->DBquad[p+k*shell->Nmodes]*qpji;
          }
          qe[m]  = qk;
          qte[m] = dqk;	  
        }
      }

      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){

          dfloat dqdr = 0, dqds = 0, dqdt = 0;
          for(int n=0;n<mesh->Nq;++n){
            dqdr += mesh->D[i*mesh->Nq+n]*qe[j*mesh->Nq + n];
            dqds += mesh->D[j*mesh->Nq+n]*qe[n*mesh->Nq + i];
          }

          int m = j*mesh->Nq + i;
          qre[m] = dqdr;
          qse[m] = dqds;
        }
      }

      interpolateQuad2D(mesh->cubInterp, qe,  mesh->Nq, cubqe,  mesh->cubNq);

      interpolateQuad2D(mesh->cubInterp, qre, mesh->Nq, cubqre, mesh->cubNq);
      interpolateQuad2D(mesh->cubInterp, qse, mesh->Nq, cubqse, mesh->cubNq);
      interpolateQuad2D(mesh->cubInterp, qte, mesh->Nq, cubqte, mesh->cubNq);

      interpolateQuad2D(mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
      interpolateQuad2D(mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);
      interpolateQuad2D(mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);

      for(int j=0;j<mesh->cubNq;++j){
        for(int i=0;i<mesh->cubNq;++i){

          int m = j*mesh->cubNq + i;

          hlong id = e*mesh->cubNp+m + k*(mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp;

          hlong gbase = e*mesh->cubNp*mesh->Nvgeo + m;
          dfloat rx = mesh->cubvgeo[gbase + RXID*mesh->cubNp];
          dfloat sx = mesh->cubvgeo[gbase + SXID*mesh->cubNp];
          dfloat ry = mesh->cubvgeo[gbase + RYID*mesh->cubNp];
          dfloat sy = mesh->cubvgeo[gbase + SYID*mesh->cubNp];
          dfloat rz = mesh->cubvgeo[gbase + RZID*mesh->cubNp];
          dfloat sz = mesh->cubvgeo[gbase + SZID*mesh->cubNp];

          dfloat rhok =  shell->Rquad[k];
          dfloat xk = cubx[m]*rhok;
          dfloat yk = cuby[m]*rhok;
          dfloat zk = cubz[m]*rhok;
          dfloat rk = sqrt(xk*xk+yk*yk+zk*zk);

          // fix chain rule here (note that coordinates scaled out by shell->Rquad[k])
          cubq[id] = cubqe[m];
          cubdqdx[id] = (rx*cubqre[m] + sx*cubqse[m])/rhok + cubqte[m]*cubx[m];
          cubdqdy[id] = (ry*cubqre[m] + sy*cubqse[m])/rhok + cubqte[m]*cuby[m];
          cubdqdz[id] = (rz*cubqre[m] + sz*cubqse[m])/rhok + cubqte[m]*cubz[m];

        }
      }
    }
  }

  free(qe); free(qre); free(qse); free(qte);
  free(cubqe); free(cubqre); free(cubqse); free(cubqte);
  free(cubx); free(cuby); free(cubz);
}
