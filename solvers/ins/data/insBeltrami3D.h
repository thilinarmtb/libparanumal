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
void insFlowField3D(bcData *bc) 
  {                             
    dfloat a = M_PI/4.f; 
    dfloat d = M_PI/2.f; 
    dfloat time = bc->time; 
    dfloat x = bc->x; 
    dfloat y = bc->y; 
    dfloat z = bc->z; 
    bc->uP = -a*exp(-p_nu*d*d*time)*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
    bc->vP = -a*exp(-p_nu*d*d*time)*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));
    bc->wP = -a*exp(-p_nu*d*d*time)*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));
    bc->pP = -0.5f*a*a*exp(-2.f*p_nu*d*d*time)*(2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+exp(2.f*a*z)+exp(2.f*a*y)+exp(2.f*a*x));
  }   
  

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
void insVelocityDirichletConditions3D(bcData *bc)
{                                       
  dfloat a = M_PI/4.f; 
    dfloat d = M_PI/2.f; 
    dfloat time = bc->time; 
    dfloat x = bc->x; 
    dfloat y = bc->y; 
    dfloat z = bc->z;                  
  if(bc->id==1){                        
    bc->uP = 0.f;                    
    bc->vP = 0.f;                    
    bc->wP = 0.f;                    
  } else if(bc->id==2){                     
    bc->uP = -a*exp(-p_nu*d*d*time)*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
    bc->vP = -a*exp(-p_nu*d*d*time)*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));
    bc->wP = -a*exp(-p_nu*d*d*time)*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));
  }else if(bc->id==4){         
    bc->uP = 0.f;       
  }else if(bc->id==5){         
    bc->vP = 0.f;       
  }else if(bc->id==6){         
    bc->wP  = 0.f;       
  }                        
}

void insPressureDirichletConditions3D(bcData *bc)
{                                   
 dfloat a = M_PI/4.f; 
    dfloat d = M_PI/2.f; 
    dfloat time = bc->time; 
    dfloat x = bc->x; 
    dfloat y = bc->y; 
    dfloat z = bc->z;                 
  if(bc->id==3){                 
    bc->pP = -0.5f*a*a*exp(-2.f*p_nu*d*d*time)*(2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+exp(2.f*a*z)+exp(2.f*a*y)+exp(2.f*a*x));
  } 
}

void insVelocityNeumannConditions3D(bcData *bc)
{                                          
  dfloat a = M_PI/4.f; 
    dfloat d = M_PI/2.f; 
    dfloat time = bc->time; 
    dfloat x = bc->x; 
    dfloat y = bc->y; 
    dfloat z = bc->z; 
  if(bc->id==3){                        
    bc->uxP = -a*(a*exp(a*x)*sin(a*y+d*z)-a*exp(a*z)*sin(a*x+d*y))*exp(-p_nu*d*d*time); 
    bc->uyP = -a*(a*exp(a*x)*cos(a*y+d*z)-d*exp(a*z)*sin(a*x+d*y))*exp(-p_nu*d*d*time); 
    bc->uzP = -a*(d*exp(a*x)*cos(a*y+d*z)+a*exp(a*z)*cos(a*x+d*y))*exp(-p_nu*d*d*time); 
    bc->vxP = -a*(d*exp(a*y)*cos(a*z+d*x)+a*exp(a*x)*cos(a*y+d*z))*exp(-p_nu*d*d*time); 
    bc->vyP = -a*(a*exp(a*y)*sin(a*z+d*x)-a*exp(a*x)*sin(a*y+d*z))*exp(-p_nu*d*d*time); 
    bc->vzP = -a*(a*exp(a*y)*cos(a*z+d*x)-d*exp(a*x)*sin(a*y+d*z))*exp(-p_nu*d*d*time); 
    bc->wxP =  a*(a*exp(a*z)*cos(a*x+d*y)-d*exp(a*y)*sin(a*z+d*x))*exp(-p_nu*d*d*time); 
    bc->wyP =  a*(d*exp(a*z)*cos(a*x+d*y)+a*exp(a*y)*cos(a*z+d*x))*exp(-p_nu*d*d*time); 
    bc->wzP =  a*(a*exp(a*z)*sin(a*x+d*y)-a*exp(a*y)*sin(a*z+d*x))*exp(-p_nu*d*d*time); 
  }                                                                                 
}

// I dont think we need this, AK. 
void insPressureNeumannConditions3D(bcData *bc)
{                                          
}


