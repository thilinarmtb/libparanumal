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
#define insFlowField2D(t,x,y,u,v,p)   \
  {                                   \
    *(u) = p_ubar;                    \
    *(v) = p_vbar;                    \
    *(p) = p_pbar;                    \
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insAdvectionDirichletConditions2D(bc, t, x, y, nx, ny, uM, vM, uB, vB) \
{                                   \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==2){                 \
    *(uB) = p_ubar;                 \
    *(vB) = p_vbar;                 \
  } else if(bc==3){                 \
    *(uB) = uM;                     \
    *(vB) = vM;                     \
  } else if(bc==4){                 \
    *(uB) = 0.f;                    \
    *(vB) = vM;                     \
  } else if(bc==5){                 \
    *(uB) = uM;                     \
    *(vB) = 0.f;                    \
  }                                 \
}


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define insVelocityDirichletConditions2D(bc, velId, t, x, y, nx, ny, uM, uB) \
{                                   \
  if(bc==1){                        \
    *(uB) = 0.f;                    \
  } else if(bc==2){                 \
    if(velId==0)                    \
      *(uB) = p_ubar;               \
    else                            \
      *(uB) = p_vbar;               \
  } else if(bc==3){                 \
    *(uB) = uM;                     \
  } else if(bc==4){                 \
    if(velId==0)                    \
    *(uB) = 0.f;                    \
    else                            \
    *(uB) = uM;                     \
  } else if(bc==5){                 \
    if(velId==0)                    \
    *(uB) = uM;                     \
    else                            \
    *(uB) = 0.f;                    \
  }                                 \
}

#define insVelocityNeumannConditions2D(bc, velId, t, x, y, nx, ny, uxM, uyM,uxB, uyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(uxB) = uxM;                          \
    *(uyB) = uyM;                          \
  } else if(bc==3){                        \
    *(uxB) = 0.f;                          \
    *(uyB) = 0.f;                          \
  } else if(bc==4){                        \
    if(velId==0){                          \
      *(uxB) = uxM;                        \
      *(uyB) = uyM;                        \
    }else{                                 \
      *(uxB) = 0.f;                        \
      *(uyB) = 0.f;                        \
    }                                      \
  } else if(bc==5){                        \
    if(velId==0){                          \
      *(uxB) = 0.f;                        \
      *(uyB) = 0.f;                        \
    }else{                                 \
      *(uxB) = uxM;                        \
      *(uyB) = uyM;                        \
    }                                      \
  }                                        \
}



#define insPressureDirichletConditions2D(bc, t, x, y, nx, ny, pM, pB) \
{                                   \
  if(bc==1 || bc==2){               \
    *(pB) = pM;                     \
  } else if(bc==3){                 \
    *(pB) = p_pbar;                 \
  } else if(bc==4){                 \
    *(pB) = pM;                     \
  } else if(bc==5){                 \
    *(pB) = pM;                     \
  }                                 \
}

#define insPressureNeumannConditions2D(bc, t, x, y, nx, ny, pxM, pyM, pxB, pyB) \
{                                          \
  if(bc==1 || bc==2){                      \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  } else if(bc==3){                        \
    *(pxB) = pxM;                          \
    *(pyB) = pyM;                          \
  } else if(bc==4){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  } else if(bc==5){                        \
    *(pxB) = 0.f;                          \
    *(pyB) = 0.f;                          \
  }                                        \
}
