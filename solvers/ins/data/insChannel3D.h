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

//    
// Initial conditions 
void cdsScalarField3D(bcData *bc) {                               
    bc->sP = 0.f;                
  }  
  
// Initial conditions 
void insFlowField3D(bcData *bc) 
  {                                  
    bc->uP = 1.f - 4.f*y*y/1.f ;    
    bc->vP = 0.f;                   
    bc->wP = 0.f;                    
    bc->pP = 0.f;                   
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
void insVelocityDirichletConditions3D(bcData *bc) 
{                                   
  if(bc->id==2){                 
    bc->uP = 1.f - 4.f*y*y/1.f ;    
    bc->vP = 0.f;                    
    bc->wP = 0.f;                    
  } 
}

void insVelocityNeumannConditions3D(bcData *bc) 
{                                          
  if(bc->id==3){                        
    bc->uxP = 0.f;                          
    bc->uyP = 0.f;                          
    bc->uzP = 0.f;                          
    bc->vxP = 0.f;                          
    bc->vyP = 0.f;                          
    bc->vzP = 0.f;                          
    bc->wxP = 0.f;                          
    bc->wyP = 0.f;                          
    bc->wzP = 0.f;                          
  }                                        
}


void insPressureDirichletConditions3D(bcData *bc) 
{                                   
  if(bc->id==3){                 
    bc->pB = 0.f;                 
  } 
}

void insPressureNeumannConditions3D(bcData *bc)
{                                          
}
