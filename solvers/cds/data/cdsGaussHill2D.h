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
void cdsFlowField2D(bcData *bc){				       
    bc->uP = -bc->y;	
    bc->vP =  bc->x;	
  }   

void cdsScalarField2D(bcData *bc){					
    dfloat mtime = bc->time + M_PI; 
    dfloat cond  = 0.001; 
    dfloat xc = 0.00f;		
    dfloat yc = 0.50f;    
    dfloat xt = xc*cos(mtime)  - yc*sin(mtime);
    dfloat yt = -xc*sin(mtime) + yc*cos(mtime);
    dfloat r2 = (bc->x-xt)*(bc->x-xt) + (bc->y-yt)*(bc->y-yt); 
    dfloat sExact =  1.f / (4.f*M_PI*cond*mtime) * exp(-r2/ (4.f*cond*mtime));
    bc->sP = sExact;				
  }   


// Boundary conditions
/* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5, z-slip 6 */
void cdsDirichletConditions2D(bcData *bc){                                   
    dfloat mtime = bc->time + M_PI/2; 
    dfloat cond  = 0.001; 
    dfloat xc = 0.00f;    
    dfloat yc = 0.50f;    
    dfloat xt = xc*cos(mtime) - yc*sin(mtime);
    dfloat yt = -xc*sin(mtime) + yc*cos(mtime);
    dfloat r2 = (bc->x-xt)*(bc->x-xt) + (bc->y-yt)*(bc->y-yt); 
    dfloat sExact =  1.f / (4.f*M_PI*cond*mtime) * exp(-r2/ (4.f*cond*mtime));  
  if(bc->id==1){                        
    bc->sP = sExact;        
  }          
}

void cdsNeumannConditions2D(bcData *bc){                                         
  if(bc->id==2){               
    bc->sxP = bc->sxM;                    
    bc->syP = bc->syM;                    
  }                                    
}
