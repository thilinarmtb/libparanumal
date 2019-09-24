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
void cdsFlowField3D(bcData *bc){				      
    bc->uP = 0.5;	
    bc->vP = 0.5; 
    bc->wP = 0.5; 
} 

// initial width = 0.05; diffusion coeff = 0.001  
void cdsScalarField3D(bcData *bc){					
    dfloat sigma_0 = 0.05; 
    dfloat sigma_d = 2.f*0.01*bc->time; 
    dfloat u = 0.5f; 
    dfloat v = 0.5f; 
    dfloat w = 0.5f; 
    dfloat xc = 0.50f;
    dfloat yc = 0.50f;
    dfloat zc = 0.50f;
    dfloat xt = xc + u*bc->time;
    dfloat yt = yc + v*bc->time;
    dfloat zt = zc + w*bc->time;
    dfloat r2 = (bc->x-xt)*(bc->x-xt) + (bc->y-yt)*(bc->y-yt) + (bc->z-zt)*(bc->z-zt);   
    dfloat sExact =  sigma_0/(sigma_0 + sigma_d) * exp(- r2 / (2.f*(sigma_0 + sigma_d)));	
    bc->sP = sExact;				
  }   


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
void cdsDirichletConditions3D(bcData *bc){                                   
    dfloat sigma_0 = 0.05; 
    dfloat sigma_d = 2.f*0.01*bc->time;
    dfloat u = 0.5f; 
    dfloat v = 0.5f; 
    dfloat w = 0.5f; 
    dfloat xc = 0.50f;
    dfloat yc = 0.50f;
    dfloat zc = 0.50f;
    dfloat xt = xc + u*bc->time;
    dfloat yt = yc + v*bc->time;
    dfloat zt = zc + w*bc->time;
    dfloat r2 = (bc->x-xt)*(bc->x-xt) + (bc->y-yt)*(bc->y-yt) + (bc->z-zt)*(bc->z-zt);
    dfloat sExact =  sigma_0/(sigma_0 + sigma_d) * exp(- r2 / (2.f*(sigma_0 + sigma_d)));
  if(bc.id==1){                        
    bc->sP = 0.f;                    
  } else if(bc.id==2){                 
    bc->sP = sExact;		    
  } else if(bc.id==3){                
    bc->sP = bc->sM;                    
  } else if(bc.id==4 || bc.id==5 || bc.id==6){  
    bc->sP = bc->sM; 
  }             
}

void cdsNeumannConditions3D(bcData *bc){                                      
 if(bc.id==1 || bc.id==2){             
    bc->sxP = bc->sxM;                      
    bc->syP = bc->syM;                      
    bc->szP = bc->szM;                      
  } else if(bc.id==3){                 
    bc->sxP = 0.f;                      
    bc->syP = 0.f;                      
    bc->szP = 0.f;                      
  } else if(bc.id==4||bc.id==5||bc.id==6){
    bc->sxP = bc->nx*bc->nx*bc->sxM;                   
    bc->syP = bc->nx*bc->nx*bc->syM;                   
    bc->szP = bc->nx*bc->nx*bc->szM;                   
  }                                       
}
