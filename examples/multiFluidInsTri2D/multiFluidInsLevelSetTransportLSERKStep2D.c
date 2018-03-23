#include "multiFluidIns2D.h"

// complete a time step using LSERK4
void multiFluidInsLevelSetTransportLSERKStep2D(multiFluidIns_t *solver, int tstep, int haloBytes,
				                          dfloat * sendBuffer, dfloat *recvBuffer, char * options){



mesh2D *mesh = solver->mesh; 

int Nfields = 1; 


// field offset at this step // index = 0 for levelSet run only, 
int offset = solver->index*(mesh->Nelements+mesh->totalHaloPairs);

// LSERK4 stages
for(int rk=0;rk<mesh->Nrk;++rk){

  // intermediate stage time
  dfloat t = tstep*solver->dt + solver->dt*mesh->rkc[rk];

  if(mesh->totalHaloPairs>0){

      solver->levelSetHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 solver->o_Phi,
                                 solver->o_phiHaloBuffer);

      // copy extracted halo to HOST
      solver->o_phiHaloBuffer.copyTo(sendBuffer);

      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*sizeof(dfloat),
                           sendBuffer,
                           recvBuffer);
    }
    
    occaTimerTic(mesh->device, "LevelSetVolumeKernel");   
    solver->levelSetTransportVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_cubDrWT,
                                          mesh->o_cubDsWT,
                                          mesh->o_cubInterpT,
                                          offset,
                                          solver->o_U,
                                          solver->o_V,
                                          solver->o_Phi,
                                          solver->o_rhsPhi);         

    occaTimerToc(mesh->device, "LevelSetVolumeKernel");   

      if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      solver->o_phiHaloBuffer.copyFrom(recvBuffer);

      solver->levelSetHaloScatterKernel(mesh->Nelements,
                                        mesh->totalHaloPairs,
                                        mesh->o_haloElementList,
                                        solver->o_Phi,
                                        solver->o_phiHaloBuffer);
      }

      occaTimerTic(mesh->device,"LevelSetSurfaceKernel");
      solver->levelSetTransportSurfaceKernel(mesh->Nelements,
                                          mesh->o_sgeo,
                                          mesh->o_intInterpT,
                                          mesh->o_intLIFTT,
                                          mesh->o_vmapM,
                                          mesh->o_vmapP,
                                          mesh->o_EToB,
                                          t,
                                          mesh->o_intx,
                                          mesh->o_inty,
                                          offset,
                                          solver->o_U,
                                          solver->o_V,
                                          solver->o_Phi,
                                          solver->o_rhsPhi);  
      
      occaTimerToc(mesh->device,"LevelSetSurfaceKernel");
 
 //  // ramp function for flow at next RK stage
 //  dfloat tupdate = tstep*mesh->dt + mesh->dt*mesh->rkc[rk+1];
 //  dfloat rampUpdate, drampdtUpdate;
 //  boltzmannRampFunction2D(tupdate, &rampUpdate, &drampdtUpdate);

 //  //UPDATE
 //  occaTimerTic(mesh->device,"UpdateKernel");

    occaTimerTic(mesh->device,"LevelSetUpdateKernel");
    solver->levelSetTransportUpdateKernel(mesh->Nelements,
                      solver->dt,
                      mesh->rka[rk],
                      mesh->rkb[rk],
                      solver->o_rhsPhi,
                      solver->o_resPhi,
                      solver->o_Phi);
    occaTimerToc(mesh->device,"LevelSetUpdateKernel");

    
 }


}
