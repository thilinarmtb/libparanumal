
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

#include <limits.h>
#include "stokes.h"

// Adapted from Method 1 (eqns 4 & 5) of
// Comput. Methods Appl. Mech. Engrg. 163 (1998) 193-204
// Projection techniques for iterative solution of Ax =b with
// successive right-hand sides
// Paul F. Fischer

void stokesFSRStart(stokes_t *stokes,
		    int &stokesNFSR, // active number of FSR rhs
		    occa::memory &b,
		    occa::memory &btilde){

  // zero accumulators
  stokes->o_fsrAlphas.copyFrom(stokes->o_fsrZeroArray);
  
  if(stokesNFSR>0){
    
    // compute Q.historyQ for all history state records
    stokes->multipleGlobalWeightedInnerProductsKernel(stokes->NtotalV,
						      stokes->ogs->o_invDegree,
						      stokes->NtotalP,
						      stokes->meshP->ogs->o_invDegree,
						      stokesNFSR,
						      b,  
						      stokes->o_fsrBtilde,
						      stokes->o_fsrAlphas);
  }
  
  // projected RHS is btilde = b - alpha_k*btilde_k
  printf("stokesNFSR=%d\n", stokesNFSR);
  stokes->fsrReconstructKernel(stokes->Ndof,
			       stokesNFSR,
			       b,
			       (dfloat) -1.0,
			       stokes->o_fsrAlphas,
			       stokes->o_fsrBtilde,
			       btilde);


  dfloat *alphas = (dfloat*) calloc(stokes->fsrNrhs, sizeof(dfloat));
  stokes->o_fsrAlphas.copyTo(alphas);
  printf("alphas = [ "); 
  for(int n=0;n<stokesNFSR;++n){
    printf("%lg ", alphas[n]);
  }
  printf("];\n");
  free(alphas);
}


void stokesFSRUpdate(stokes_t *stokes,
		     int &stokesNFSR, // active number of FSR rhs
		     dfloat lambda,
		     occa::memory &xtilde,
		     occa::memory &x,
		     occa::memory &btilde){

  // x = xtilde + alpha_k*xtilde_k
  stokes->fsrReconstructKernel(stokes->Ndof,
			       stokesNFSR,
			       xtilde,
			       (dfloat) +1.0,
			       stokes->o_fsrAlphas,
			       stokes->o_fsrXtilde,
			       x);
  
  if(stokesNFSR==stokes->fsrNrhs){
    
    stokesOperator(stokes, lambda, x, btilde);
    
    dfloat normbtilde2 = 0;
    stokesVecNorm2(stokes, btilde, &normbtilde2);
    dfloat invNormbtilde = 1./sqrt(normbtilde2);
    
    stokes->vecScaledAddKernel(stokes->Ndof, invNormbtilde, btilde, 0, stokes->o_fsrBtilde);
    stokes->vecScaledAddKernel(stokes->Ndof, invNormbtilde,      x, 0, stokes->o_fsrXtilde);

    stokesNFSR = 1;
  }
  else{

    occa::memory bhat = btilde;
    
    stokesOperator(stokes, lambda, xtilde, bhat);
    
    // zero accumulators
    stokes->o_fsrAlphas.copyFrom(stokes->o_fsrZeroArray);
    
    // compute Axtilde.Btilde for all history state records
    stokes->multipleGlobalWeightedInnerProductsKernel(stokes->NtotalV,
						      stokes->ogs->o_invDegree,
						      stokes->NtotalP,
						      stokes->meshP->ogs->o_invDegree,
						      stokesNFSR,
						      bhat,  
						      stokes->o_fsrBtilde,
						      stokes->o_fsrAlphas);


    // bhat  - alpha_k*Btilde_k
    stokes->fsrReconstructKernel(stokes->Ndof,
				 stokesNFSR,
				 bhat,
				 (dfloat) -1.0,
				 stokes->o_fsrAlphas,
				 stokes->o_fsrBtilde,
				 btilde);

    
    // x  - alpha_k*Xtilde_k
    stokes->fsrReconstructKernel(stokes->Ndof,
				 stokesNFSR,
				 xtilde,
				 (dfloat) -1.0,
				 stokes->o_fsrAlphas,
				 stokes->o_fsrXtilde,
				 xtilde);

    // ||Ax  - alpha_k*Btilde_k||

    dfloat normbtilde2 = 0;
    stokesVecNorm2(stokes, btilde, &normbtilde2);
    dfloat invNormbtilde = 1./sqrt(normbtilde2);
    
    stokes->vecScaledAddKernel(stokes->Ndof, invNormbtilde, btilde, 0,
			       stokes->o_fsrBtilde + stokesNFSR*stokes->Ndof*sizeof(dfloat));
    
    stokes->vecScaledAddKernel(stokes->Ndof, invNormbtilde, xtilde, 0,
			       stokes->o_fsrXtilde + stokesNFSR*stokes->Ndof*sizeof(dfloat));
    
    ++stokesNFSR;

  }
}
