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
#include "omp.h"
#include <unistd.h>

static dfloat asbfManufacturedForcingFunction(asbf_t *asbf, dfloat x, dfloat y, dfloat z);

asbf_t *asbfSetup(mesh_t *mesh, dfloat lambda, occa::properties kernelInfo, setupAide options){

  asbf_t *asbf = (asbf_t*) calloc(1, sizeof(asbf_t));
  asbf->mesh = mesh;
  asbf->options = options;

  options.getArgs("MESH DIMENSION", asbf->dim);
  options.getArgs("ELEMENT TYPE", asbf->elementType);

  mesh->Nfields = 1;

  if(asbf->dim==3){
    if(asbf->elementType == QUADRILATERALS){
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
    }else if(asbf->elementType == TRIANGLES){
      meshOccaSetupTri3D(mesh, options, kernelInfo);
    }else{
      meshOccaSetup3D(mesh, options, kernelInfo);
    }
  }
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  asbfSolveSetup(asbf, lambda, kernelInfo);

  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat xbase = mesh->x[e*mesh->Np+n];
      dfloat ybase = mesh->y[e*mesh->Np+n];
      dfloat zbase = mesh->z[e*mesh->Np+n];

      dfloat J;

      if(asbf->elementType==QUADRILATERALS)
        J = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JID) + n];
      else
        J = mesh->vgeo[e*mesh->Nvgeo + JID];

      for(int g=0;g<asbf->Nquad;++g){

        dfloat Rg = asbf->Rquad[g];

        // stretch coordinates
        dfloat xg = Rg*xbase;
        dfloat yg = Rg*ybase;
        dfloat zg = Rg*zbase;

        // evaluate rhs at asbf quadrature for each surface node
        asbf->f[g] = asbfManufacturedForcingFunction(asbf, xg, yg, zg);
      }

      // integrate f against asbf modes
      for(int m=0;m<asbf->Nmodes;++m){
        dfloat fhatm = 0;
        for(int i=0;i<asbf->Nquad;++i){
          fhatm += asbf->Bquad[m + i*asbf->Nmodes]*asbf->Wquad[i]*asbf->f[i];
        }

        // scale by surface weight
        asbf->r3D[e*mesh->Np + n + m*asbf->Ntotal] = J*fhatm;
      }
    }
  }


  return asbf;
}

/******************************************************************************/

// NB:  This must match asbfManufacturedSolution() in asbfErrorHex3D.c
static dfloat asbfManufacturedForcingFunction(asbf_t *asbf, dfloat x, dfloat y, dfloat z)
{
  dfloat r, theta, phi;
  dfloat f;

  r     = sqrt(x*x + y*y + z*z);
  theta = atan2(y, x);
  phi   = acos(z/r);

  dfloat k1 = 6.283185307179586;
  dfloat k2 = 18.849555921538759;
  dfloat k3 = 25.132741228718345;
  f = (k1 + asbf->lambda/k1)*sin(k1*r)/r
          + (k2 + asbf->lambda/k2)*sin(k2*r)/r
          + (k3 + asbf->lambda/k3)*sin(k3*r)/r;

  return f;
}
