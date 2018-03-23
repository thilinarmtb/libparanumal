#include "multiFluidIns2D.h"

void multiFluidInsLevelSetRun2D(multiFluidIns_t *levelSet, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = levelSet->mesh; 
  
  // MPI send buffer
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  int haloBytes;

  haloBytes = mesh->totalHaloPairs*mesh->Np*levelSet->Nfields*sizeof(dfloat);
  
  if (haloBytes) {
    occa::memory o_sendBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBufferPinned = mesh->device.mappedAlloc(haloBytes, NULL);
    sendBuffer = (dfloat*) o_sendBufferPinned.getMappedPointer();
    recvBuffer = (dfloat*) o_recvBufferPinned.getMappedPointer();
  }
  
// printf("N: %d Nsteps: %d dt: %.5e \n", mesh->N, mesh->NtimeSteps, mesh->dt);

  for(int tstep=0;tstep<levelSet->NtimeSteps;++tstep){

   

    // if(strstr(options, "LSERK")){
    //     // occaTimerTic(mesh->device, "LSERK");  
        multiFluidInsLevelSetTransportLSERKStep2D(levelSet, tstep, haloBytes, sendBuffer, recvBuffer, options);
    //     // occaTimerToc(mesh->device, "LSERK");  
    //   }


    if(strstr(options, "REPORT")){
      if((tstep%levelSet->errorStep)==0){
        multiFluidInsReport2D(levelSet, tstep, options);
      }
    }


  }



printf("writing Final data\n");  
multiFluidInsReport2D(levelSet, levelSet->NtimeSteps, options);

// occa::printTimer();
}


