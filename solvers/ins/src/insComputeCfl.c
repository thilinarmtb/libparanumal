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

#include "ins.h"

void insGetDh(ins_t *ins); 

dfloat insComputeCfl(ins_t *ins, dfloat time, int tstep){

  mesh_t *mesh = ins->mesh; 
  // create dH i.e. nodal spacing in reference coordinates
  if(!ins->cflComputed)
    insGetDh(ins);
  // Compute cfl factors i.e. dt* U / h 
  ins->cflKernel(mesh->Nelements,
                 ins->dt, 
                 mesh->o_vgeo,
                 ins->o_idH,
                 ins->fieldOffset,
                 ins->o_U,
                 ins->o_rhsU);  
  
  // find the local maximum of CFL number
  const dlong Ntotal = mesh->Np*mesh->Nelements; 
  dfloat *tmp        = (dfloat *) calloc(ins->Nblock, sizeof(dfloat));
  occa::memory o_tmp = mesh->device.malloc(ins->Nblock*sizeof(dfloat), tmp);
  ins->maxKernel(Ntotal, ins->o_rhsU, o_tmp);
  o_tmp.copyTo(tmp);
  
  // finish reduction
  dfloat cfl = 0.f; 
  for(dlong n=0;n<ins->Nblock;++n){
    cfl  = mymax(cfl, tmp[n]);
  }

  dfloat gcfl = 0.f;
  MPI_Allreduce(&cfl, &gcfl, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
  
  free(tmp);
  return gcfl; 
}


void insGetDh(ins_t *ins){

  mesh_t *mesh = ins->mesh; 
  
  dfloat *dH; 
  if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA){
    dH = (dfloat*) calloc((mesh->N+1),sizeof(dfloat));

    for(int n=0; n<(mesh->N+1); n++){
      if(n==0)
        dH[n] = mesh->gllz[n+1] - mesh->gllz[n];
      else if(n==mesh->N)
        dH[n] = mesh->gllz[n] - mesh->gllz[n-1];
      else
        dH[n] = 0.5*( mesh->gllz[n+1] - mesh->gllz[n-1]); 
    }
    // // invert it
    for(int n=0; n< (mesh->N+1); n++)
      dH[n] = 1.0/dH[n]; 

    // Move data to Device
    ins->o_idH = mesh->device.malloc((mesh->N+1)*sizeof(dfloat), dH); 
    free(dH); 
  }else{
    printf("cfl for tri and tet is not included yet\n");
    exit(0);
  }

  ins->cflComputed = 1; 

}
