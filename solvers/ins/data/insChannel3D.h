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
#define insFlowField3D(t,x,y,z, u,v,w,p) \
  {                                   \
    *(u) = 1.f - 4.f*y*y/1.f ;     \
    *(v) = 0.f;                    \
    *(w) = 0.f;                    \
    *(p) = 0.f;                    \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define insVelocityDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, uM, vM, wM, uB, vB, wB) \
{                                   \
  if(bc==2){                 \
    *(uB) = 1.f - 4.f*y*y/1.f ;     \
    *(vB) = 0.f;                    \
    *(wB) = 0.f;                    \
  } \
}

#define insVelocityNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, uxM, uyM, uzM, vxM, vyM, vzM, wxM, wyM, wzM, uxB, uyB, uzB, vxB, vyB, vzB, wxB, wyB, wzB) \
{                                          \
  if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
    *(uzB) = 0.f;                          \
    *(vxB) = 0.f;                          \
    *(vyB) = 0.f;                          \
    *(vzB) = 0.f;                          \
    *(wxB) = 0.f;                          \
    *(wyB) = 0.f;                          \
    *(wzB) = 0.f;                          \
  }                                        \
}


#define insPressureDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, pM, pB) \
{                                   \
  if(bc==3){                 \
    *(pB) = 0.f;                 \
  } \
}

#define insPressureNeumannConditions3D(bc, t, x, y, z, nx, ny, nz, pxM, pyM, pzM, pxB, pyB, pzB) \
{                                          \
}
