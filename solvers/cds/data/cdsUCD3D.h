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
void cdsFlowField3D(bcData *bc) {				       
    bc->uP = 0.8f;
    bc->vP = 0.8f;
    bc->wP = 0.8f;
  }   

void cdsScalarField3D(bcData *bc){					
    dfloat ax = 0.01f;
    dfloat ay = 0.01f;
    dfloat az = 0.01f;
    dfloat cx = 0.8f;
    dfloat cy = 0.8f;
    dfloat cz = 0.8f;
    dfloat coef   = 1.f/pow((4.f*bc->time +1.f),1.5f);			
    dfloat xterm  = pow((bc->x - cx*bc->time - 0.5f),2.f) / ( ax*(4.f*bc->time+1.f) );
    dfloat yterm  = pow((bc->y - cy*bc->time - 0.5f),2.f) / ( ay*(4.f*bc->time+1.f) );	
    dfloat zterm  = pow((bc->z - cz*bc->time - 0.5f),2.f) / ( az*(4.f*bc->time+1.f) );	
    dfloat sExact = coef*exp( -xterm -yterm -zterm);
    bc->sP = sExact;
  }   


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
void cdsDirichletConditions3D(bcData *bc){                                   
    dfloat ax = 0.01f;			
    dfloat ay = 0.01f;			
    dfloat az = 0.01f;			
    					
    dfloat cx = 0.8f;
    dfloat cy = 0.8f;
    dfloat cz = 0.8f;
    						
    dfloat coef  = 1.f/pow((4.f*bc->time +1.f),1.5f);
    dfloat xterm = pow((bc->x - cx*bc->time- 0.5f),2.f) / ( ax*(4.f*bc->time+1.f) );	
    dfloat yterm = pow((bc->y - cy*bc->time- 0.5f),2.f) / ( ay*(4.f*bc->time+1.f) );	
    dfloat zterm = pow((bc->z - cz*bc->time- 0.5f),2.f) / ( az*(4.f*bc->time+1.f) );	
    									
    dfloat sExact = coef*exp(-xterm -yterm - zterm);
  if(bc->id==1){  
    bc->sP = sExact;
  }       
}

void cdsNeumannConditions3D(bcData *bc){          
  if(bc->id==2){                  
    bc->sxP = bc->sxM;                      
    bc->syP = bc->syM;                      
    bc->szP = bc->szM;                      
  }                             
}
