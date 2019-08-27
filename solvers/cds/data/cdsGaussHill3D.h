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


// Initial conditions 
#define cdsFlowField3D(t,x,y,z,u,v,w) \
  {				       \
    *(u) = 0.5;		       \
    *(v) = 0.5;           \
    *(w) = 0.5;           \
  } 

// initial width = 0.05; diffusion coeff = 0.001  
#define cdsScalarField3D(t,x,y,z,s)	\
  {					\
    dfloat sigma_0 = 0.05; \
    dfloat sigma_d = 2.f*0.01*t; \
    dfloat u = 0.5f; \
    dfloat v = 0.5f; \
    dfloat w = 0.5f; \
    dfloat xc = 0.50f;			\
    dfloat yc = 0.50f;     \
    dfloat zc = 0.50f;     \
    dfloat xt = xc + u*t;      \
    dfloat yt = yc + v*t;      \
    dfloat zt = zc + w*t;      \
    dfloat r2 = (x-xt)*(x-xt) + (y-yt)*(y-yt) + (z-zt)*(z-zt);      \
    dfloat sExact =  sigma_0/(sigma_0 + sigma_d) * exp(- r2 / (2.f*(sigma_0 + sigma_d)));	\
    *(s) = sExact;				\
  }   


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define cdsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, sM,sB) \
{                                   \
    dfloat sigma_0 = 0.05; \
    dfloat sigma_d = 2.f*0.01*t; \
    dfloat u = 0.5f; \
    dfloat v = 0.5f; \
    dfloat w = 0.5f; \
    dfloat xc = 0.50f;      \
    dfloat yc = 0.50f;     \
    dfloat zc = 0.50f;     \
    dfloat xt = xc + u*t;      \
    dfloat yt = yc + v*t;      \
    dfloat zt = zc + w*t;      \
    dfloat r2 = (x-xt)*(x-xt) + (y-yt)*(y-yt) + (z-zt)*(z-zt);      \
    dfloat sExact =  sigma_0/(sigma_0 + sigma_d) * exp(- r2 / (2.f*(sigma_0 + sigma_d))); \
  if(bc==1){                        \
    *(sB) = 0.f;                    \
  } else if(bc==2){                 \
    *(sB) = sExact;		    \
  } else if(bc==3){                 \
    *(sB) = sM;                     \
  } else if(bc==4||bc==5||bc==6){   \
    *(sB) = sM; \
  }                                 \
}

#define cdsNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, sxM, syM,szM, sxB, syB, szB) \
{                                          \
 if(bc==1 || bc==2){                      \
    *(sxB) = sxM;                          \
    *(syB) = syM;                          \
    *(szB) = szM;                          \
  } else if(bc==3){                        \
    *(sxB) = 0.f;                          \
    *(syB) = 0.f;                          \
    *(szB) = 0.f;                          \
  } else if(bc==4||bc==5||bc==6){          \
    *(sxB) = nx*nx*sxM;                    \
    *(syB) = nx*nx*syM;                    \
    *(szB) = nx*nx*szM;                    \
  }                                        \                    \
}
