
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

void stokesFSRStart(stokes_t *stokes,
		    int &stokesNFSR, // active number of FSR rhs
		    const stokesVec_t &b,
		    stokesVec_t &u){
  
  // zero accumulators
  stokes->o_fsrAlphas.copyFrom(stokes->o_fsrZeroArray);
  
  // compute Q.historyQ for all history state records
  stokes->multipleGlobalWeightedInnerProductsKernel(stokes->NtotalV,
						    stokes->ogs->o_invDegree,
						    stokes->NtotalP,
						    stokes->meshP->ogs->o_invDegree,
						    stokesNFSR,
						    b.o_v,  
						    stokes->o_fsrHistoryQ,
						    stokes->o_fsrAlphas);

  // new intial guess is in u
  stokes->fsrStartKernel(stokes->Ndof,
			 stokesNFSR,
			 stokes->o_fsrAlphas,
			 stokes->o_fsrHistoryQ,
			 u.o_v);
  
}


void stokesFSRUpdate(stokes_t *stokes,
		     int &stokesNFSR, // active number of FSR rhs
		     dfloat lambda,
		     stokesVec_t &v,
		     stokesVec_t &Av
		     ){

  // apply stokes Operator
  stokesOperator(stokes, lambda, v, Av);

  dfloat vdotAv = 0;
  stokesVecInnerProduct(stokes, v, Av, &vdotAv);
  
  if(stokesNFSR==stokes->fsrNrhs-1){
    
    v.o_v.copyTo(stokes->o_fsrHistoryQ, stokes->Ndof*sizeof(dfloat), 0); // 0 offset
    
    dfloat fac = 1./vdotAv;
    
    stokesVecScale(stokes, v, fac);
  }
  else{
    // zero accumulators
    stokes->o_fsrAlphas.copyFrom(stokes->o_fsrZeroArray);
    
    // compute Q.historyQ for all history state records
    stokes->multipleGlobalWeightedInnerProductsKernel(stokes->NtotalV,
						      stokes->ogs->o_invDegree,
						      stokes->NtotalP,
						      stokes->meshP->ogs->o_invDegree,
						      stokesNFSR,
						      Av.o_v,  
						      stokes->o_fsrHistoryQ,
						      stokes->o_fsrAlphas);

    
    dfloat *alphas = (dfloat*) calloc(stokesNFSR, sizeof(dfloat));
    stokes->o_fsrAlphas.copyTo(alphas, stokesNFSR*sizeof(dfloat));

    dfloat normAlphas2 = 0;
    for(int fld=0;fld<stokesNFSR;++fld){
      normAlphas2 += pow(alphas[fld],2);
    }      

    dfloat invNormAlphas = 1./sqrt(normAlphas2);

    stokes->fsrUpdateKernel(stokes->Ndof,
			    stokesNFSR,
			    invNormAlphas,
			    stokes->o_fsrAlphas,
			    v.o_v,
			    stokes->o_fsrHistoryQ);
    
    ++stokesNFSR;
  }
}
