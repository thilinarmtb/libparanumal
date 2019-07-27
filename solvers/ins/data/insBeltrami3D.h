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
#define insFlowField3D(bc) \
  {                                   \
    dfloat a = M_PI/4.f; \
    dfloat d = M_PI/2.f; \
    *(bc.uP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.x)*sin(a*bc.y+d*bc.z)+exp(a*bc.z)*cos(a*bc.x+d*bc.y));\
    *(bc.vP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.y)*sin(a*bc.z+d*bc.x)+exp(a*bc.x)*cos(a*bc.y+d*bc.z));\
    *(bc.wP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.z)*sin(a*bc.x+d*bc.y)+exp(a*bc.y)*cos(a*bc.z+d*bc.x));\
    *(bc.pP) = -0.5f*a*a*exp(-2.f*p_nu*d*d*bc.time)*(2.f*exp(a*(bc.z+bc.y))*cos(a*bc.z+d*bc.x)*sin(a*bc.x+d*bc.y)+2.f*exp(a*(bc.z+bc.x))*cos(a*bc.x+d*bc.y)*sin(a*bc.y+d*bc.z)+2.f*exp(a*(bc.y+bc.x))*cos(a*bc.y+d*bc.z)*sin(a*bc.z+d*bc.x)+exp(2.f*a*bc.z)+exp(2.f*a*bc.y)+exp(2.f*a*bc.x));\
  }   
  

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, bc.z-slip 6 */
#define insVelocityDirichletConditions3D(bc) \
{                                       \
  dfloat a = M_PI/4.f;                  \
  dfloat d = M_PI/2.f;                  \
  if(bc.id==1){                        \
    *(bc.uP) = 0.f;                    \
    *(bc.vP) = 0.f;                    \
    *(bc.wP) = 0.f;                    \
  } else if(bc.id==2){                     \
    *(bc.uP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.x)*sin(a*bc.y+d*bc.z)+exp(a*bc.z)*cos(a*bc.x+d*bc.y));\
    *(bc.vP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.y)*sin(a*bc.z+d*bc.x)+exp(a*bc.x)*cos(a*bc.y+d*bc.z));\
    *(bc.wP) = -a*exp(-p_nu*d*d*bc.time)*(exp(a*bc.z)*sin(a*bc.x+d*bc.y)+exp(a*bc.y)*cos(a*bc.z+d*bc.x));\
  }else if(bc.id==4){         \
    *(bc.uP) = 0.f;       \
  }else if(bc.id==5){         \
    *(bc.vP) = 0.f;       \
  }else if(bc.id==6){         \
    *(bc.wP)  = 0.f;       \
  }                        \
}

#define insPressureDirichletConditions3D(bc) \
{                                   \
  dfloat a = M_PI/4.f;                  \
  dfloat d = M_PI/2.f;                  \
  if(bc.id==3){                 \
    *(bc.pP) = -0.5f*a*a*exp(-2.f*p_nu*d*d*bc.time)*(2.f*exp(a*(bc.z+bc.y))*cos(a*bc.z+d*bc.x)*sin(a*bc.x+d*bc.y)+2.f*exp(a*(bc.z+bc.x))*cos(a*bc.x+d*bc.y)*sin(a*bc.y+d*bc.z)+2.f*exp(a*(bc.y+bc.x))*cos(a*bc.y+d*bc.z)*sin(a*bc.z+d*bc.x)+exp(2.f*a*bc.z)+exp(2.f*a*bc.y)+exp(2.f*a*bc.x));\
  } \
}

#define insVelocityNeumannConditions3D(bc) \
{                                          \
  dfloat a = M_PI/4.f; \
  dfloat d = M_PI/2.f; \
  if(bc.id==3){                        \
    *(bc.uxP) = -a*(a*exp(a*bc.x)*sin(a*bc.y+d*bc.z)-a*exp(a*bc.z)*sin(a*bc.x+d*bc.y))*exp(-p_nu*d*d*bc.time); \
    *(bc.uyP) = -a*(a*exp(a*bc.x)*cos(a*bc.y+d*bc.z)-d*exp(a*bc.z)*sin(a*bc.x+d*bc.y))*exp(-p_nu*d*d*bc.time); \
    *(bc.uzP) = -a*(d*exp(a*bc.x)*cos(a*bc.y+d*bc.z)+a*exp(a*bc.z)*cos(a*bc.x+d*bc.y))*exp(-p_nu*d*d*bc.time); \
    *(bc.vxP) = -a*(d*exp(a*bc.y)*cos(a*bc.z+d*bc.x)+a*exp(a*bc.x)*cos(a*bc.y+d*bc.z))*exp(-p_nu*d*d*bc.time); \
    *(bc.vyP) = -a*(a*exp(a*bc.y)*sin(a*bc.z+d*bc.x)-a*exp(a*bc.x)*sin(a*bc.y+d*bc.z))*exp(-p_nu*d*d*bc.time); \
    *(bc.vzP) = -a*(a*exp(a*bc.y)*cos(a*bc.z+d*bc.x)-d*exp(a*bc.x)*sin(a*bc.y+d*bc.z))*exp(-p_nu*d*d*bc.time); \
    *(bc.wxP) =  a*(a*exp(a*bc.z)*cos(a*bc.x+d*bc.y)-d*exp(a*bc.y)*sin(a*bc.z+d*bc.x))*exp(-p_nu*d*d*bc.time); \
    *(bc.wyP) =  a*(d*exp(a*bc.z)*cos(a*bc.x+d*bc.y)+a*exp(a*bc.y)*cos(a*bc.z+d*bc.x))*exp(-p_nu*d*d*bc.time); \
    *(bc.wzP) =  a*(a*exp(a*bc.z)*sin(a*bc.x+d*bc.y)-a*exp(a*bc.y)*sin(a*bc.z+d*bc.x))*exp(-p_nu*d*d*bc.time); \
  }                                                                                 \
}

// I dont think we need this, AK. 
#define insPressureNeumannConditions3D(bc) \
{                                          \
}


