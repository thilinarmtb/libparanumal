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

#include "elliptic.h"
#include "ogsInterface.h"

#include "omp.h"


//host gs
void ellipticSerialOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision){
  
  mesh_t *mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);
  
  //  printf("serialOperator: gathers = %d, reductions = %d, cts = %d, serial = %d, ipdg = %d\n",
  //	 enableGatherScatters, enableReductions, continuous, serial, ipdg);
  
  dfloat *sendBuffer = elliptic->sendBuffer;
  dfloat *recvBuffer = elliptic->recvBuffer;
  dfloat *gradSendBuffer = elliptic->gradSendBuffer;
  dfloat *gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;

  if(continuous && serial){

    ellipticSerialAxHexKernel3D(mesh->Nq,  mesh->Nelements, mesh->o_ggeo, 
				mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda, o_q, o_Aq, elliptic->o_ggeoNoJW);
    
    ogs_t *ellipticOgs = elliptic->ogs;
    
    ogsHostGatherScatter(o_Aq.ptr(), dfloatString, "add", ellipticOgs->hostGsh);

#if USE_NULL_BOOST==1
    if(elliptic->allNeumann) { // inspect this later
      
      dlong Nblock = elliptic->Nblock;
      dfloat *tmp = elliptic->tmp;
      occa::memory &o_tmp = elliptic->o_tmp;

      elliptic->innerProductKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }
#endif
    
    //post-mask
    if (elliptic->Nmasked) 
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  } else if(ipdg){
    printf("WARNING: DEBUGGING C0\n");
    MPI_Finalize();
    exit(-1);
  } 
  
}

void ellipticOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq,
	const char *precision){

  mesh_t *mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);

  //  printf("generalOperator: gathers = %d, reductions = %d, cts = %d, serial = %d, ipdg = %d\n",
  //	 enableGatherScatters, enableReductions, continuous, serial, ipdg);
  
  dfloat *sendBuffer = elliptic->sendBuffer;
  dfloat *recvBuffer = elliptic->recvBuffer;
  dfloat *gradSendBuffer = elliptic->gradSendBuffer;
  dfloat *gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = elliptic->Nblock;
  dfloat *tmp = elliptic->tmp;
  occa::memory &o_tmp = elliptic->o_tmp;

  if(continuous){

    // TW: turned off for debugging
    if(serial && elliptic->elementType==HEXAHEDRA){
      ellipticSerialOperator(elliptic, lambda, o_q, o_Aq, precision);
      return;
    }

    ogs_t *ogs = elliptic->ogs;

    int mapType = (elliptic->elementType==HEXAHEDRA &&
                   options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;
    
    int integrationType = (elliptic->elementType==HEXAHEDRA &&
			   options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;
    
    occa::kernel &partialAxKernel =
	(strstr(precision, "float")) ? elliptic->partialFloatAxKernel : elliptic->partialAxKernel;
    
    if(mesh->NglobalGatherElements) {
      
      if(integrationType==0) { // GLL or non-hex
	if(mapType==0){
	  partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
			  mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, lambda,
			  o_q, o_Aq);
	}
	else
	  partialAxKernel(mesh->NglobalGatherElements, mesh->o_globalGatherElementList,
			  elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices,
			  mesh->o_MM, lambda, o_q, o_Aq);
      }
      else{
	elliptic->partialCubatureAxKernel(mesh->NglobalGatherElements,
					  mesh->o_globalGatherElementList,
					  mesh->o_cubggeo,
					  mesh->o_cubD,
					  mesh->o_cubInterpT,
					  lambda, o_q, o_Aq);
      }
    }
    
    if(enableGatherScatters)
      ogsGatherScatterStart(o_Aq, ogsDfloat, ogsAdd, ogs);
    
    if(mesh->NlocalGatherElements){
      if(integrationType==0) { // GLL or non-hex
	if(mapType==0){
	  partialAxKernel(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
			  mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM,
		      	  lambda, o_q, o_Aq);
	}
	else
	  partialAxKernel(mesh->NlocalGatherElements, mesh->o_localGatherElementList,
			  elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices,
			  mesh->o_MM, lambda, o_q, o_Aq);
      }
      else{
	elliptic->partialCubatureAxKernel(mesh->NlocalGatherElements,
					  mesh->o_localGatherElementList,
					  mesh->o_cubggeo,
					  mesh->o_cubD,
					  mesh->o_cubInterpT,
					  lambda,
					  o_q,
					  o_Aq);
      }
    }

    // finalize gather using local and global contributions
    if(enableGatherScatters==1)
      ogsGatherScatterFinish(o_Aq, ogsDfloat, ogsAdd, ogs);

#if USE_NULL_BOOST==1
    if(elliptic->allNeumann) {

      elliptic->innerProductKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;

      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
    }
#endif
    
    //post-mask
    if (elliptic->Nmasked) 
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  } else if(ipdg){
    dlong offset = 0;
    dfloat alpha = 0., alphaG =0.;
    dlong Nblock = elliptic->Nblock;
    dfloat *tmp = elliptic->tmp;
    occa::memory &o_tmp = elliptic->o_tmp;

    ellipticStartHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    if(options.compareArgs("BASIS", "NODAL")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_Dmatrices,
				      o_q,
				      elliptic->o_grad);
    } else if(options.compareArgs("BASIS", "BERN")) {
      elliptic->partialGradientKernel(mesh->Nelements,
				      offset,
				      mesh->o_vgeo,
				      mesh->o_D1ids,
				      mesh->o_D2ids,
				      mesh->o_D3ids,
				      mesh->o_Dvals,
				      o_q,
				      elliptic->o_grad);
    }

    ellipticInterimHaloExchange(elliptic, o_q, mesh->Np, sendBuffer, recvBuffer);

    //Start the rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
#if USE_NULL_BOOST==1
    if(elliptic->allNeumann)
      mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
#endif
    
    if(mesh->NinternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialIpdgKernel(mesh->NinternalElements,
				    mesh->o_internalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }

#if USE_NULL_BOOST==1
    if(elliptic->allNeumann) {
      o_tmp.copyTo(tmp);

      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
      alphaG *= elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
    }
#endif
    
    ellipticEndHaloExchange(elliptic, o_q, mesh->Np, recvBuffer);

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_Dmatrices,
					o_q,
					elliptic->o_grad);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialGradientKernel(mesh->totalHaloPairs,
					offset,
					mesh->o_vgeo,
					mesh->o_D1ids,
					mesh->o_D2ids,
					mesh->o_D3ids,
					mesh->o_Dvals,
					o_q,
					elliptic->o_grad);
      }
    }

    if(mesh->NnotInternalElements) {
      if(options.compareArgs("BASIS", "NODAL")) {
        elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_Dmatrices,
				    mesh->o_LIFTT,
				    mesh->o_MM,
				    elliptic->o_grad,
				    o_Aq);
      } else if(options.compareArgs("BASIS", "BERN")) {
        elliptic->partialIpdgKernel(mesh->NnotInternalElements,
				    mesh->o_notInternalElementIds,
				    mesh->o_vmapM,
				    mesh->o_vmapP,
				    lambda,
				    elliptic->tau,
				    mesh->o_vgeo,
				    mesh->o_sgeo,
				    elliptic->o_EToB,
				    mesh->o_D1ids,
				    mesh->o_D2ids,
				    mesh->o_D3ids,
				    mesh->o_Dvals,
				    mesh->o_L0vals,
				    mesh->o_ELids,
				    mesh->o_ELvals,
				    mesh->o_BBMM,
				    elliptic->o_grad,
				    o_Aq);
      }
    }

#if USE_NULL_BOOST==1
    if(elliptic->allNeumann)
      mesh->addScalarKernel(mesh->Nelements*mesh->Np, alphaG, o_Aq);
#endif
  } 

}

