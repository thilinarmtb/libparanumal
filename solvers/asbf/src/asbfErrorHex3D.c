#include "asbf.h"

void asbfErrorHex3D(mesh_t *mesh, dfloat *q){
  
  dfloat normErrorH1 = 0;
  dfloat normErrorL2 = 0;
  
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
	  
	  dfloat K = 6.283185307179586;
	  //dfloat k = 18.849555921538759;
	  //dfloat k = 25.132741228718345;
	  dfloat qE = sin(K*rg)/(K*rg);
	  dfloat dqEdrho = cos(K*rg)/rg - sin(K*rg)/(K*rg*rg); // rho is radial

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

  normErrorH1 = sqrt(normErrorH1);
  normErrorL2 = sqrt(normErrorL2);
  
  printf("%g, %g %%%% abs H1 error norm, L2 error norm\n", normErrorH1, normErrorL2);
  
}
