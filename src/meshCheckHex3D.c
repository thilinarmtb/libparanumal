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

#include "mesh3D.h"

void meshCheckHex3D(mesh_t *mesh){

  dfloat vtol = 1e-9;

  for(hlong e=0;e<mesh->Nelements;++e){
    for(hlong n=0;n<mesh->Np;++n){
      hlong gid = e*mesh->Nvgeo*mesh->Np + n + JID*mesh->Np;
      if(mesh->vgeo[gid]<vtol){
	printf("small Jacobian: element %d, node %d, J=%lg\n",
	       e, n, mesh->vgeo[gid]);
      }
    }
  }

  dfloat stol = 1e-6;
	  
  // check surface geometric factors
  for(hlong eM=0;eM<mesh->Nelements;++eM){
    for(int fM=0;fM<mesh->Nfaces;++fM){
      for(int nM=0;nM<mesh->Nfp;++nM){

	hlong idM = eM*mesh->Nfaces*mesh->Nfp + fM*mesh->Nfp + nM;

	hlong vidM = mesh->vmapM[idM];
	hlong sidM  = mesh->Nsgeo*idM;

	dfloat nxM = mesh->sgeo[sidM+NXID];
	dfloat nyM = mesh->sgeo[sidM+NYID];
	dfloat nzM = mesh->sgeo[sidM+NZID];
	dfloat sJM = mesh->sgeo[sidM+SJID];

	dfloat xM = mesh->x[vidM];
	dfloat yM = mesh->y[vidM];
	dfloat zM = mesh->z[vidM];
	
	if(sJM<stol){
	  printf("small sJ: (%d,%d,%d, %g) ", eM, fM, nM, sJM);
	}
	
	hlong idP = mesh->mapP[idM];
	
	if(idP>=0){
	  hlong eP = idP/(mesh->Nfaces*mesh->Nfp);
	  hlong fP = (idP%(mesh->Nfaces*mesh->Nfp))/mesh->Nfp;
	  hlong nP = (idP%mesh->Nfp);
	  hlong vidP = mesh->vmapP[idM];
	  
	  hlong sidP  = mesh->Nsgeo*idP;
	  
	  dfloat nxP = mesh->sgeo[sidP+NXID];
	  dfloat nyP = mesh->sgeo[sidP+NYID];
	  dfloat nzP = mesh->sgeo[sidP+NZID];
	  dfloat sJP = mesh->sgeo[sidP+SJID];

	  dfloat xP = mesh->x[vidP];
	  dfloat yP = mesh->y[vidP];
	  dfloat zP = mesh->z[vidP];

	  dfloat coordDiscrepancy =
	    fabs(xM-xP) + fabs(yM-yP) + fabs(zM-zP);

	  if(coordDiscrepancy > stol)
	    printf("node coordinate mismatch: (%d,%d,%d, %g,%g,%g) "
		   " => (%d,%d,%d, %g,%g,%g)\n",
		   eM, fM, nM, xM, yM, zM,
		   eP, fP, nP, xP, yP, zP);
	  
	  if(vidM!=vidP){
	    dfloat normalDiscrepancy =
	      fabs(nxM+nxP) + fabs(nyM+nyP) + fabs(nzM+nzP);
	    
	    dfloat sJDiscrepancy = fabs(sJM-sJP);
	    
	    if(normalDiscrepancy > stol)
	      printf("normal mismatch: (%d,%d,%d, %g,%g,%g) "
		     " => (%d,%d,%d, %g,%g,%g)\n",
		     eM, fM, nM, nxM, nyM, nzM,
		     eP, fP, nP, nxP, nyP, nzP);
	    
	    if(sJDiscrepancy > stol)
	      printf("surface Jacobian mismatch: (%d,%d,%d, %g) "
		     " => (%d,%d,%d, %g)\n",
		     eM, fM, nM, sJM,
		     eP, fP, nP, sJP);
	  }
	}
      }
    }
  }
}
