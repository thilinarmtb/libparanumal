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

extern void interpolateQuad2D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);
extern void interpolateHex3D(dfloat *I, dfloat *x, int N, dfloat *Ix, int M);

static void shellErrorHex3DASBF(shell_t *shell, dfloat *q3D, dfloat *normErrorH1, dfloat *normErrorL2);
static void shellErrorHex3DSEM(shell_t *shell, dfloat *q3D, dfloat *normErrorH1, dfloat *normErrorL2);

static void shellManufacturedSolution(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz);
static void shellManufacturedSolutionDD(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz);
static void shellManufacturedSolutionNN(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz);

void shellErrorHex3D(shell_t *shell, dfloat *q3D, dfloat *normErrorH1, dfloat *normErrorL2)
{
  setupAide options = shell->options;

  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    shellErrorHex3DASBF(shell, q3D, normErrorH1, normErrorL2);
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    shellErrorHex3DSEM(shell, q3D, normErrorH1, normErrorL2);
  } else {
    printf("ERROR:  Invalid value \"%s\" for SHELL SOLVER.\n",
           options.getArgs("SHELL SOLVER").c_str());
    exit(-1);
  }
}

static void shellErrorHex3DASBF(shell_t *shell, dfloat *q3D, dfloat *normErrorH1, dfloat *normErrorL2)
{
  mesh_t *mesh = shell->mesh;

  *normErrorH1 = 0;
  *normErrorL2 = 0;

  hlong cubNtotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp*shell->Nquad;
  dfloat *cubq    = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdx = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdy = (dfloat*) calloc(cubNtotal, sizeof(dfloat));
  dfloat *cubdqdz = (dfloat*) calloc(cubNtotal, sizeof(dfloat));

  dfloat* cubx = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cuby = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubz = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  shellCubatureGradient(shell, q3D, cubq, cubdqdx, cubdqdy, cubdqdz);

  for(int e = 0; e < mesh->Nelements; e++){
    interpolateQuad2D(mesh->cubInterp, mesh->x+e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
    interpolateQuad2D(mesh->cubInterp, mesh->y+e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);
    interpolateQuad2D(mesh->cubInterp, mesh->z+e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);

    for(int k=0;k<shell->Nquad;++k){
      for(int j=0;j<mesh->cubNq;++j){
        for(int i=0;i<mesh->cubNq;++i){
          int m = j*mesh->cubNq + i;

          hlong id = e*mesh->cubNp+m + k*(mesh->Nelements+mesh->totalHaloPairs)*mesh->cubNp;

          dfloat qg    = cubq[id];
          dfloat dqdxg = cubdqdx[id];
          dfloat dqdyg = cubdqdy[id];
          dfloat dqdzg = cubdqdz[id];

          hlong gbase = e*mesh->cubNp*mesh->Nvgeo + m;
          dfloat JW = mesh->cubvgeo[gbase + JWID*mesh->cubNp]*shell->Wquad[k];

          dfloat rhog = shell->Rquad[k];

          dfloat xg = cubx[m]*rhog;
          dfloat yg = cuby[m]*rhog;
          dfloat zg = cubz[m]*rhog;

          dfloat qE, dqEdx, dqEdy, dqEdz;
          shellManufacturedSolution(shell, xg, yg, zg, &qE, &dqEdx, &dqEdy, &dqEdz);

          dfloat localH1 = pow(dqEdx-dqdxg, 2) + pow(dqEdy-dqdyg, 2) + pow(dqEdz-dqdzg, 2) + pow(qE-qg, 2);
          dfloat localL2 = pow(qE-qg, 2);

          *normErrorH1 += JW*localH1;
          *normErrorL2 += JW*localL2;
        }
      }
    }
  }

  free(cubq);
  free(cubdqdx); free(cubdqdy);   free(cubdqdz);

  *normErrorH1 = sqrt(*normErrorH1);
  *normErrorL2 = sqrt(*normErrorL2);
}

static void shellErrorHex3DSEM(shell_t *shell, dfloat *q3D, dfloat *normErrorH1, dfloat *normErrorL2)
{
  elliptic_t *elliptic = shell->elliptic;
  mesh_t *mesh         = elliptic->mesh;

  // TODO:  Move the gradient calculation to shellCubatureGradient().
  dfloat *qr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *qt = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *cubx = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cuby = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubz = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));

  dfloat *cubq = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqr = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqs = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *cubqt = (dfloat*)calloc(mesh->cubNp, sizeof(dfloat));

  *normErrorH1 = 0.0;
  *normErrorL2 = 0.0;
  for (int e = 0; e < mesh->Nelements; e++) {
    for (int i = 0; i < mesh->Nq; i++) {
      for (int j = 0; j < mesh->Nq; j++) {
        for (int k = 0; k < mesh->Nq; k++) {
          int ind = i*mesh->Nq*mesh->Nq + j*mesh->Nq + k;

          dfloat dqdr = 0.0;
          dfloat dqds = 0.0;
          dfloat dqdt = 0.0;

          for (int n = 0; n < mesh->Nq; n++) {
            dqdr += mesh->D[k*mesh->Nq + n]*mesh->q[e*mesh->Np + i*mesh->Nq*mesh->Nq + j*mesh->Nq + n];
            dqds += mesh->D[j*mesh->Nq + n]*mesh->q[e*mesh->Np + i*mesh->Nq*mesh->Nq + n*mesh->Nq + k];
            dqdt += mesh->D[i*mesh->Nq + n]*mesh->q[e*mesh->Np + n*mesh->Nq*mesh->Nq + j*mesh->Nq + k];
          }

          qr[ind] = dqdr;
          qs[ind] = dqds;
          qt[ind] = dqdt;
        }
      }
    }

    interpolateHex3D(mesh->cubInterp, mesh->q + e*mesh->Np, mesh->Nq, cubq, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->x + e*mesh->Np, mesh->Nq, cubx, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, mesh->y + e*mesh->Np, mesh->Nq, cuby, mesh->cubNq);

    interpolateHex3D(mesh->cubInterp, mesh->z + e*mesh->Np, mesh->Nq, cubz, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, qr, mesh->Nq, cubqr, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, qs, mesh->Nq, cubqs, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, qt, mesh->Nq, cubqt, mesh->cubNq);

    for (int i = 0; i < mesh->cubNq; i++) {
      for (int j = 0; j < mesh->cubNq; j++) {
        for (int k = 0; k < mesh->cubNq; k++) {
          hlong ind = i*mesh->cubNq*mesh->cubNq + j*mesh->cubNq + k;
          dfloat x = cubx[ind];
          dfloat y = cuby[ind];
          dfloat z = cubz[ind];

          dfloat qE, dqEdx, dqEdy, dqEdz;
          shellManufacturedSolution(shell, x, y, z, &qE, &dqEdx, &dqEdy, &dqEdz);

          hlong gind = e*mesh->cubNp*mesh->Nvgeo + i*mesh->cubNq*mesh->cubNq + j*mesh->cubNq + k;

          dfloat drdx = mesh->cubvgeo[gind + RXID*mesh->cubNp];
          dfloat dsdx = mesh->cubvgeo[gind + SXID*mesh->cubNp];
          dfloat dtdx = mesh->cubvgeo[gind + TXID*mesh->cubNp];

          dfloat drdy = mesh->cubvgeo[gind + RYID*mesh->cubNp];
          dfloat dsdy = mesh->cubvgeo[gind + SYID*mesh->cubNp];
          dfloat dtdy = mesh->cubvgeo[gind + TYID*mesh->cubNp];

          dfloat drdz = mesh->cubvgeo[gind + RZID*mesh->cubNp];
          dfloat dsdz = mesh->cubvgeo[gind + SZID*mesh->cubNp];
          dfloat dtdz = mesh->cubvgeo[gind + TZID*mesh->cubNp];

          dfloat JW = mesh->cubvgeo[gind + JWID*mesh->cubNp];

          dfloat dqdx = drdx*cubqr[ind] + dsdx*cubqs[ind] + dtdx*cubqt[ind];
          dfloat dqdy = drdy*cubqr[ind] + dsdy*cubqs[ind] + dtdy*cubqt[ind];
          dfloat dqdz = drdz*cubqr[ind] + dsdz*cubqs[ind] + dtdz*cubqt[ind];

          dfloat errH1local = pow(dqEdx - dqdx, 2) + pow(dqEdy - dqdy, 2) + pow(dqEdz - dqdz, 2) + pow(qE - cubq[ind], 2);
          dfloat errL2local = pow(qE - cubq[ind], 2);

          *normErrorH1 += JW*errH1local;
          *normErrorL2 += JW*errL2local;
        }
      }
    }
  }

  *normErrorH1 = sqrt(*normErrorH1);
  *normErrorL2 = sqrt(*normErrorL2);
}

/******************************************************************************/

static void shellManufacturedSolution(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz)
{
  if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
    return shellManufacturedSolutionDD(shell, x, y, z, q, dqdx, dqdy, dqdz);
  } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
    return shellManufacturedSolutionNN(shell, x, y, z, q, dqdx, dqdy, dqdz);
  } else {
    printf("ERROR:  No manufactured solution for type-(%d, %d) BCs.\n",
           shell->innerBC, shell->outerBC);
    exit(-1);
  }

  return;
}

// Appropriate BCs:  Dirichlet-Dirichlet
static void shellManufacturedSolutionDD(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz)
{
  dfloat r, theta, phi;
  dfloat drdx, drdy, drdz;
  dfloat dthetadx, dthetady, dthetadz;
  dfloat dphidx, dphidy, dphidz;

  r     = sqrt(x*x + y*y + z*z);
  theta = atan2(y, x);
  phi   = acos(z/r);

  drdx = cos(theta)*sin(phi);
  drdy = sin(theta)*sin(phi);
  drdz = cos(phi);

  dthetadx = -sin(theta)/(r*sin(phi));
  dthetady = cos(theta)/(r*sin(phi));
  dthetadz = 0;

  dphidx = cos(theta)*cos(phi)/r;
  dphidy = sin(theta)*cos(phi)/r;
  dphidz = -sin(phi)/r;

#if 0
  dfloat k1 = 6.283185307179586;
  dfloat k2 = 18.849555921538759;
  dfloat k3 = 25.132741228718345;
  *q = sin(k1*r)/(k1*r) + sin(k2*r)/(k2*r) + sin(k3*r)/(k3*r);
  dfloat dqdr = cos(k1*r)/r - sin(k1*r)/(k1*r*r)
    + cos(k2*r)/r - sin(k2*r)/(k2*r*r)
    + cos(k3*r)/r - sin(k3*r)/(k3*r*r);
  *dqdx = dqdr*x/r;
  *dqdy = dqdr*y/r;
  *dqdz = dqdr*z/r;
#elif 1
  dfloat r2mR2      = pow(r, 2.0) - pow(shell->R, 2.0);
  dfloat r2m1       = pow(r, 2.0) - 1.0;
  dfloat twor2mR2m1 = 2.0*pow(r, 2.0) - pow(shell->R, 2.0) - 1.0;

  *q = sin(x)*cos(y)*exp(z)*r2mR2*r2m1;
  *dqdx = cos(y)*exp(z)*(2.0*x*sin(x)*twor2mR2m1 + cos(x)*r2mR2*r2m1); 
  *dqdy = sin(x)*exp(z)*(2.0*y*cos(y)*twor2mR2m1 - sin(y)*r2mR2*r2m1);
  *dqdz = sin(x)*cos(y)*exp(z)*(2.0*z*twor2mR2m1 + r2mR2*r2m1);
#endif

  return;
}

// Appropriate BCs:  Neumann-Neumann
static void shellManufacturedSolutionNN(shell_t *shell,
                                     dfloat x, dfloat y, dfloat z, dfloat* q,
                                     dfloat* dqdx, dfloat* dqdy, dfloat* dqdz)
{
  dfloat r, theta, phi;
  dfloat drdx, drdy, drdz;
  dfloat dthetadx, dthetady, dthetadz;
  dfloat dphidx, dphidy, dphidz;

  r     = sqrt(x*x + y*y + z*z);
  theta = atan2(y, x);
  phi   = acos(z/r);

  drdx = cos(theta)*sin(phi);
  drdy = sin(theta)*sin(phi);
  drdz = cos(phi);

  dthetadx = -sin(theta)/(r*sin(phi));
  dthetady = cos(theta)/(r*sin(phi));
  dthetadz = 0;

  dphidx = cos(theta)*cos(phi)/r;
  dphidy = sin(theta)*cos(phi)/r;
  dphidz = -sin(phi)/r;

  // NB:  This solution may only work for shell->R not too much bigger than 1.
  dfloat dqdr, dqdtheta, dqdphi;
  dfloat dqdthetadthetadx, dqdthetadthetady;
  dfloat p, dpdr, s, dsdphi;

#if 0
  p = r*(pow(r, 2.0)/3.0 - ((1.0 + shell->R)/2.0)*r + shell->R);
  dpdr = (1.0 - r)*(shell->R - r);
#else
  p = cos(M_PI*(r - 1.0)/(shell->R - 1.0));
  dpdr = -M_PI*sin(M_PI*(r - 1.0)/(shell->R - 1.0))/(shell->R - 1.0);
#endif

  s = pow(phi, 3.0)*pow(phi - M_PI, 3.0);
  dsdphi = 3.0*pow(phi, 2.0)*pow(phi - M_PI, 2.0)*(2.0*phi - M_PI);

  *q = sin(theta)*s*p;

  dqdr = sin(theta)*s*dpdr;
  dqdtheta = cos(theta)*s*p;
  dqdphi = sin(theta)*dsdphi*p;

  if ((fabs(phi) < 1.0e-13) || (fabs(phi - M_PI) < 1.0e-13)) {
    // Near poles, the Laplacian is approximately zero.
    dqdthetadthetadx = 0;
    dqdthetadthetady = 0;
  } else {
    dqdthetadthetadx = dqdtheta*dthetadx;
    dqdthetadthetady = dqdtheta*dthetady;
  }

  *dqdx = dqdr*drdx + dqdthetadthetadx + dqdphi*dphidx;
  *dqdy = dqdr*drdy + dqdthetadthetady + dqdphi*dphidy;
  *dqdz = dqdr*drdz + dqdphi*dphidz;

  return;
}
