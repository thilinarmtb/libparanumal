#include "asbf.h"

void asbfErrorHex3D(asbf_t *asbf, dfloat *q3D){

  mesh_t *mesh = asbf->mesh;

  dfloat normErrorH1 = 0;
  dfloat normErrorL2 = 0;

#if 0
  for (int e = 0; e < mesh->Nelements; e++){
    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){

          hlong id = e*mesh->Np +  k*mesh->Nq*mesh->Nq + j*mesh->Nq + i;

          dfloat qg = mesh->q[id];

          dfloat xg = mesh->x[id];
          dfloat yg = mesh->y[id];
          dfloat zg = mesh->z[id];
          dfloat rg = sqrt(xg*xg + yg*yg + zg*zg);

          dfloat K1 = 6.283185307179586;
          dfloat K2 = 18.849555921538759;
          dfloat K3 = 25.132741228718345;
          dfloat qE = sin(K1*rg)/(K1*rg) + sin(K2*rg)/(K2*rg) + sin(K3*rg)/(K3*rg);
          dfloat dqEdrho = cos(K1*rg)/rg - sin(K1*rg)/(K1*rg*rg) // rho is radial
            + cos(K2*rg)/rg - sin(K2*rg)/(K2*rg*rg)
            + cos(K3*rg)/rg - sin(K3*rg)/(K3*rg*rg);

          dfloat dqEdx = dqEdrho*xg/rg;
          dfloat dqEdy = dqEdrho*yg/rg;
          dfloat dqEdz = dqEdrho*zg/rg;

          dfloat dqdr = 0, dqds = 0, dqdt = 0;
          for(int n=0;n<mesh->Nq;++n){
            dqdr += mesh->D[i*mesh->Nq+n]*q[e*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + n];
            dqds += mesh->D[j*mesh->Nq+n]*q[e*mesh->Np + k*mesh->Nq*mesh->Nq + n*mesh->Nq + i];
            dqdt += mesh->D[k*mesh->Nq+n]*q[e*mesh->Np + n*mesh->Nq*mesh->Nq + j*mesh->Nq + i];
          }

          hlong gid = e*mesh->Np*mesh->Nvgeo + k*mesh->Nq*mesh->Nq + j*mesh->Nq + i;

          dfloat drdx = mesh->vgeo[gid + RXID*mesh->Np];
          dfloat dsdx = mesh->vgeo[gid + SXID*mesh->Np];
          dfloat dtdx = mesh->vgeo[gid + TXID*mesh->Np];

          dfloat drdy = mesh->vgeo[gid + RYID*mesh->Np];
          dfloat dsdy = mesh->vgeo[gid + SYID*mesh->Np];
          dfloat dtdy = mesh->vgeo[gid + TYID*mesh->Np];

          dfloat drdz = mesh->vgeo[gid + RZID*mesh->Np];
          dfloat dsdz = mesh->vgeo[gid + SZID*mesh->Np];
          dfloat dtdz = mesh->vgeo[gid + TZID*mesh->Np];

          dfloat JW = mesh->vgeo[gid + JWID*mesh->Np];

          dfloat dqdx = drdx*dqdr + dsdx*dqds + dtdx*dqdt;
          dfloat dqdy = drdy*dqdr + dsdy*dqds + dtdy*dqdt;
          dfloat dqdz = drdz*dqdr + dsdz*dqds + dtdz*dqdt;

          dfloat localH1 = pow(dqEdx-dqdx, 2) + pow(dqEdy-dqdy, 2) + pow(dqEdz-dqdz, 2) + pow(qE-qg, 2);
          dfloat localL2 = pow(qE-qg, 2);

          normErrorH1 += JW*localH1;
          normErrorL2 += JW*localL2;

        }
      }
    }
  }
#endif

#if 0
  dfloat *cubq  = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqr = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqs = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqt = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubx  = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cuby  = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubz  = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  dfloat *qr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qt = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  for(int e = 0; e < mesh->Nelements; e++){

    // differentiate
    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){

          int id = k*mesh->Nq*mesh->Nq + j*mesh->Nq + i;

          dfloat dqdr = 0, dqds = 0, dqdt = 0;
          for(int n=0;n<mesh->Nq;++n){
            dqdr += mesh->D[i*mesh->Nq+n]*q[e*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + n];
            dqds += mesh->D[j*mesh->Nq+n]*q[e*mesh->Np + k*mesh->Nq*mesh->Nq + n*mesh->Nq + i];
            dqdt += mesh->D[k*mesh->Nq+n]*q[e*mesh->Np + n*mesh->Nq*mesh->Nq + j*mesh->Nq + i];
          }

          qr[id] = dqdr;
          qs[id] = dqds;
          qt[id] = dqdt;
        }
      }
    }

    void interpolateHex3D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);
    interpolateHex3D(mesh->cubInterp, mesh->q+e*mesh->Np, mesh->Nq, cubq, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);

    interpolateHex3D(mesh->cubInterp, qr, mesh->Nq, cubqr, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, qs, mesh->Nq, cubqs, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, qt, mesh->Nq, cubqt, mesh->cubNq);

    for(int k=0;k<mesh->cubNq;++k){
      for(int j=0;j<mesh->cubNq;++j){
        for(int i=0;i<mesh->cubNq;++i){

          hlong id = k*mesh->cubNq*mesh->cubNq + j*mesh->cubNq + i;

          dfloat cubqg = cubq[id];
          dfloat cubqrg = cubqr[id];
          dfloat cubqsg = cubqs[id];
          dfloat cubqtg = cubqt[id];

          dfloat xg = cubx[id];
          dfloat yg = cuby[id];
          dfloat zg = cubz[id];
          dfloat rg = sqrt(xg*xg + yg*yg + zg*zg);

          dfloat K1 = 6.283185307179586;
          dfloat K2 = 18.849555921538759;
          dfloat K3 = 25.132741228718345;
          dfloat qE = sin(K1*rg)/(K1*rg) + sin(K2*rg)/(K2*rg) + sin(K3*rg)/(K3*rg);
          dfloat dqEdrho = cos(K1*rg)/rg - sin(K1*rg)/(K1*rg*rg) // rho is radial
            + cos(K2*rg)/rg - sin(K2*rg)/(K2*rg*rg)
            + cos(K3*rg)/rg - sin(K3*rg)/(K3*rg*rg);

          dfloat dqEdx = dqEdrho*xg/rg;
          dfloat dqEdy = dqEdrho*yg/rg;
          dfloat dqEdz = dqEdrho*zg/rg;

          hlong gid = e*mesh->cubNp*mesh->Nvgeo + k*mesh->cubNq*mesh->cubNq + j*mesh->cubNq + i;

          dfloat drdx = mesh->cubvgeo[gid + RXID*mesh->cubNp];
          dfloat dsdx = mesh->cubvgeo[gid + SXID*mesh->cubNp];
          dfloat dtdx = mesh->cubvgeo[gid + TXID*mesh->cubNp];

          dfloat drdy = mesh->cubvgeo[gid + RYID*mesh->cubNp];
          dfloat dsdy = mesh->cubvgeo[gid + SYID*mesh->cubNp];
          dfloat dtdy = mesh->cubvgeo[gid + TYID*mesh->cubNp];

          dfloat drdz = mesh->cubvgeo[gid + RZID*mesh->cubNp];
          dfloat dsdz = mesh->cubvgeo[gid + SZID*mesh->cubNp];
          dfloat dtdz = mesh->cubvgeo[gid + TZID*mesh->cubNp];

          dfloat JW = mesh->cubvgeo[gid + JWID*mesh->cubNp];

          dfloat dqdx = drdx*cubqrg + dsdx*cubqsg + dtdx*cubqtg;
          dfloat dqdy = drdy*cubqrg + dsdy*cubqsg + dtdy*cubqtg;
          dfloat dqdz = drdz*cubqrg + dsdz*cubqsg + dtdz*cubqtg;

          dfloat localH1 = pow(dqEdx-dqdx, 2) + pow(dqEdy-dqdy, 2) + pow(dqEdz-dqdz, 2) + pow(qE-cubqg, 2);
          dfloat localL2 = pow(qE-cubqg, 2);

          normErrorH1 += JW*localH1;
          normErrorL2 += JW*localL2;

        }
      }
    }
  }

  free(qr); free(qs); free(qt);
  free(cubq);
  free(cubqr); free(cubqs);   free(cubqt);
  free(cubx); free(cuby);   free(cubz);

#endif

#if 1

  hlong cubNtotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp*asbf->asbfNquad;
  dfloat *cubq    = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdx = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdy = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdz = (dfloat*) calloc(cubNtotal, sizeof(dfloat));

  dfloat* cubx = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cuby = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubz = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  asbfCubatureGradient(asbf, q3D, cubq, cubdqdx, cubdqdy, cubdqdz);

  for(int e = 0; e < mesh->Nelements; e++){

    interpolateQuad2D(mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
    interpolateQuad2D(mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);
    interpolateQuad2D(mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);

    for(int k=0;k<asbf->asbfNquad;++k){
      for(int j=0;j<mesh->cubNq;++j){
        for(int i=0;i<mesh->cubNq;++i){

          int m = j*mesh->cubNq + i;

          hlong id = e*mesh->cubNp+m + k*(mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp;

          dfloat qg    = cubq[id];
          dfloat dqdxg = cubdqdx[id];
          dfloat dqdyg = cubdqdy[id];
          dfloat dqdzg = cubdqdz[id];

          hlong gbase = e*mesh->cubNp*mesh->Nvgeo + m;
          dfloat JW = mesh->cubvgeo[gbase + JWID*mesh->cubNp]*asbf->asbfWquad[k];

          dfloat rhog = asbf->asbfRquad[k];

          dfloat xg = cubx[m]*rhog;
          dfloat yg = cuby[m]*rhog;
          dfloat zg = cubz[m]*rhog;
          dfloat rg = sqrt(xg*xg + yg*yg + zg*zg);

          dfloat K1 = 6.283185307179586;
          dfloat K2 = 18.849555921538759;
          dfloat K3 = 25.132741228718345;
          dfloat qE = sin(K1*rg)/(K1*rg) + sin(K2*rg)/(K2*rg) + sin(K3*rg)/(K3*rg);
          dfloat dqEdrho = cos(K1*rg)/rg - sin(K1*rg)/(K1*rg*rg) // rho is radial
            + cos(K2*rg)/rg - sin(K2*rg)/(K2*rg*rg)
            + cos(K3*rg)/rg - sin(K3*rg)/(K3*rg*rg);

          dfloat dqEdx = dqEdrho*xg/rg;
          dfloat dqEdy = dqEdrho*yg/rg;
          dfloat dqEdz = dqEdrho*zg/rg;

          dfloat localH1 = pow(dqEdx-dqdxg, 2) + pow(dqEdy-dqdyg, 2) + pow(dqEdz-dqdzg, 2) + pow(qE-qg, 2);
          dfloat localL2 = pow(qE-qg, 2);

          //printf("dqdx = %g,%g - dqdy = %g,%g - dqdz = %g,%g - q = %g,%g - JW = %g\n",
          //       dqEdx, dqdxg, dqEdy, dqdyg, dqEdz, dqdzg, qE, qg, JW);

          normErrorH1 += JW*localH1;
          normErrorL2 += JW*localL2;

        }
      }
    }
  }

  free(cubq);
  free(cubdqdx); free(cubdqdy);   free(cubdqdz);

#endif



  normErrorH1 = sqrt(normErrorH1);
  normErrorL2 = sqrt(normErrorL2);

  printf("%g, %g %%%% abs H1 error norm, L2 error norm\n", normErrorH1, normErrorL2);

}
