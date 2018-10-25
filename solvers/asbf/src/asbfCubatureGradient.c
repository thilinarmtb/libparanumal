#include "asbf.h"

void interpolateQuad2D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M){

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

void asbfCubatureGradient(asbf_t *asbf, dfloat *q3D,
    dfloat *cubq, dfloat *cubdqdx, dfloat *cubdqdy, dfloat *cubdqdz){

  mesh_t *mesh = asbf->mesh;

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
    for(int k=0;k<asbf->asbfNquad;++k){

      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){

          int m = j*mesh->Nq + i;
          dfloat qk = 0, dqk = 0;
          for(int p=0;p<asbf->asbfNmodes;++p){
            dfloat qpji = asbf->q3D[(e*mesh->Np+m)+p*asbf->Ntotal];
            qk  +=  asbf->asbfBquad[p+k*asbf->asbfNmodes]*qpji;
            dqk += asbf->asbfDBquad[p+k*asbf->asbfNmodes]*qpji;
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

          dfloat rhok =  asbf->asbfRquad[k];
          dfloat xk = cubx[m]*rhok;
          dfloat yk = cuby[m]*rhok;
          dfloat zk = cubz[m]*rhok;
          dfloat rk = sqrt(xk*xk+yk*yk+zk*zk);

          // fix chain rule here (note that coordinates scaled out by asbf->asbfRquad[k])
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
