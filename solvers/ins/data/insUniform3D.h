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
#define insFlowField3D(bcData *bc) 
  {                                   
    bc->uP = p_ubar;                    
    bc->vP = p_vbar;                    
    bc->wP = p_wbar;                    
    bc->pP = p_pbar;                    
  }   

// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
#define insVelocityDirichletConditions3D(bcData *bc)
{                                   
 if(bc->id==2){                 
    bc->uP = p_ubar;                 
    bc->vP = p_vbar;                 
    bc->wP = p_wbar;                 
  } else if(bc->id==3){                 
    bc->uP = bc->uM;                     
    bc->vP = bc->vM;                     
    bc->wP = bc->wM;                     
  }                              
}

#define insVelocityNeumannConditions3D(bcData *bc)
{                                          
  if(bc->id==1 || bc->id==2){                      
    bc->uxP = bc->uxM;                          
    bc->uyP = bc->uyM;                          
    bc->uzP = bc->uzM;                          
    bc->vxP = bc->vxM;                          
    bc->vyP = bc->vyM;                          
    bc->vzP = bc->vzM;                          
    bc->wxP = bc->wxM;                          
    bc->wyP = bc->wyM;                          
    bc->wzP = bc->wzM;                          
  } else if(bc->id==3){                        
    bc->uxP = 0.f;                          
    bc->uyP = 0.f;                          
    bc->uzP = 0.f;                          
    bc->vxP = 0.f;                          
    bc->vyP = 0.f;                          
    bc->vzP = 0.f;                          
    bc->wxP = 0.f;                          
    bc->wyP = 0.f;                          
    bc->wzP = 0.f;                          
  } else if(bc->id==4||bc->id==5||bc->id==6){          
    bc->uxP = bc->nx*bc->nx*bc->uxM;                    
    bc->uyP = bc->nx*bc->nx*bc->uyM;                    
    bc->uzP = bc->nx*bc->nx*bc->uzM;                    
    bc->vxP = bc->ny*bc->ny*bc->vxM;                    
    bc->vyP = bc->ny*bc->ny*bc->vyM;                    
    bc->vzP = bc->ny*bc->ny*bc->vzM;                    
    bc->wxP = bc->nz*bc->nz*bc->wxM;                    
    bc->wyP = bc->nz*bc->nz*bc->wyM;                    
    bc->wzP = bc->nz*bc->nz*bc->wzM;                    
  }                                        
}


#define insPressureDirichletConditions3D(bcData *bc) 
{                                   
  if(bc->id==3){                 
    bc->pP = p_pbar;                 
  }                            
}

#define insPressureNeumannConditions3D(bcData *bc)
{                                          
  if(bc->id==3){                        
    bc->pxP = bc->pxM;                          
    bc->pyP = bc->pyM;                          
    bc->pzP = bc->pzM;                          
  }                                      
}
