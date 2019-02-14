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

#include "../stokes.h"

int main(int argc, char **argv)
{
  stokes_t         *stokes;
  occa::properties kernelInfoV, kernelInfoP;
  dfloat           *p, *pRaised;
  occa::memory     o_p, o_pRaised, o_interpRaise;

  MPI_Init(&argc, &argv);

  setupAide options("testRaisePressureHex3D.rc");

  stokes = stokesSetup(kernelInfoP, kernelInfoV, options);

  /*
  for (int e = 0; e < stokes->meshP->Nelements; e++) {
    printf("Element %d:  ", e);
    for (int i = 0; i < stokes->meshP->Nverts; i++) {
      printf("%d ", stokes->meshP->EToV[e*stokes->meshP->Nverts + i]);
    }
    printf("\n");
  }

  for (int e = 0; e < stokes->meshP->Nelements; e++) {
    printf("Element %d:\n", e);
    for (int i = 0; i < stokes->meshP->Np; i++) {
      dfloat r, s, t;
      dfloat x, y, z;
      dfloat rx, ry, rz;
      dfloat sx, sy, sz;
      dfloat tx, ty, tz;
      dfloat J;

      r = stokes->meshP->r[i];
      s = stokes->meshP->s[i];
      t = stokes->meshP->t[i];

      x = stokes->meshP->x[i];
      y = stokes->meshP->y[i];
      z = stokes->meshP->z[i];

      rx = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*RXID + i];
      ry = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*RYID + i];
      rz = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*RZID + i];
      sx = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*SXID + i];
      sy = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*SYID + i];
      sz = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*SZID + i];
      tx = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*TXID + i];
      ty = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*TYID + i];
      tz = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*TZID + i];
      J  = stokes->meshP->vgeo[e*stokes->meshP->Np*stokes->meshP->Nvgeo + stokes->meshP->Np*JID  + i];

      printf("  Point %d:  (r, s, t) --> (x, y, z) : (% .15e, % .15e, % .15e) --> (% .15e, % .15e, % .15e)\n", i, r, s, t, x, y, z);
      printf("  RX = % .15e,  RY = % .15e,  RZ = % .15e\n", rx, ry, rz);
      printf("  SX = % .15e,  SY = % .15e,  SZ = % .15e\n", sx, sy, sz);
      printf("  TX = % .15e,  TY = % .15e,  TZ = % .15e\n", tx, ty, tz);
      printf("  J  = % .15e\n", J);

    }
  }
  */

  p = (dfloat*)calloc(stokes->NtotalP, sizeof(dfloat));

  printf("Reference:\n");
  for (int i = 0; i < stokes->meshP->Np; i++) {
    dfloat r, s, t;

    r = stokes->meshP->r[i];
    s = stokes->meshP->s[i];
    t = stokes->meshP->t[i];

    printf("(% .15e, % .15e, % .15e)\n", r, s, t);
  }

  printf("Before:\n");
  for (int i = 0; i < stokes->NtotalP; i++) {
    dfloat x, y, z;

    x = stokes->meshP->x[i];
    y = stokes->meshP->y[i];
    z = stokes->meshP->z[i];
    p[i] = (double)i;

    printf("(% .15e, % .15e, % .15e) --> % .15e\n", x, y, z, p[i]);
  }
  o_p = stokes->meshV->device.malloc(stokes->NtotalP*sizeof(dfloat), p);

  pRaised = (dfloat*)calloc(stokes->NtotalV, sizeof(dfloat));
  o_pRaised = stokes->meshV->device.malloc(stokes->NtotalV*sizeof(dfloat), pRaised);

  o_interpRaise = stokes->meshV->device.malloc(stokes->meshP->Nq*stokes->meshV->Nq*sizeof(dfloat), stokes->meshP->interpRaise);

  stokes->raisePressureKernel(stokes->meshV->Nelements,
                              o_interpRaise,
                              o_p,
                              o_pRaised);

  o_pRaised.copyTo(pRaised);

  printf("After:\n");
  for (int i = 0; i < stokes->NtotalV; i++) {
    dfloat x, y, z;

    x = stokes->meshV->x[i];
    y = stokes->meshV->y[i];
    z = stokes->meshV->z[i];

    printf("(% .15e, % .15e, % .15e) --> % .15e\n", x, y, z, pRaised[i]);
  }

  free(p);
  free(pRaised);

  MPI_Finalize();

  return 0;
}
