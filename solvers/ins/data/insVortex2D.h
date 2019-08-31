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
void insFlowField2D(bcData *bc) 
  {                             
    bc->uP = -sin(2.f*M_PI*bc->y)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
    bc->vP =  sin(2.f*M_PI*bc->x)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
    bc->pP = -cos(2.f*M_PI*bc->y)*cos(2.f*M_PI*bc->x)*exp(-p_nu*8.f*M_PI*M_PI*bc->time); 
  }   
// Default u+ = u-, modify if it is different  
void insVelocityDirichletConditions2D(bcData *bc)
{                                                       
  if(bc->id==2){                                            
    bc->uP = -sin(2.f*M_PI*bc->y)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
    bc->vP =  sin(2.f*M_PI*bc->x)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
  }
}

void insVelocityNeumannConditions2D(bcData *bc)
{                                         
  if(bc->id==3){                               
    bc->uxP = 0.f;                          
    bc->uyP =-2.f*M_PI*cos(2.f*M_PI*bc->y)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
    bc->vxP = 2.f*M_PI*cos(2.f*M_PI*bc->x)*exp(-p_nu*4.f*M_PI*M_PI*bc->time);
    bc->vyP = 0.f;                          
  }                                        
}

void insPressureDirichletConditions2D(bcData *bc)
{                                  
  if(bc->id==3){                        
    bc->pP = -cos(2.f*M_PI*bc->y)*cos(2.f*M_PI*bc->x)*exp(-p_nu*8.f*M_PI*M_PI*bc->time);
   } 
}

void insPressureNeumannConditions2D(bcData *bc)
{                                          
  if(bc->id==3){                               
  }                                        
}