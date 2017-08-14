#include "ellipticHex3D.h"

occa::kernel saferBuildKernelFromSource(occa::device &device, 
					const char *fname, const char *kname, occa::kernelInfo &kernelInfo){
  
  // should really use root to build and non-root to load
  return device.buildKernelFromSource(fname, kname, kernelInfo);
  
}



// specialized version for geometric factors at Gauss (not GLL) nodes
dfloat *ellipticGeometricFactorsHex3D(mesh3D *mesh){

  /* number of second order geometric factors */
  iint NgjGeo = 7;
  iint gjNq = mesh->gjNq;
  iint gjNp = gjNq*gjNq*gjNq;
  dfloat *gjGeo = (dfloat*) calloc(mesh->Nelements*NgjGeo*gjNp, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    int id = e*mesh->Nverts;
    
    dfloat *xe = mesh->EX + id;
    dfloat *ye = mesh->EY + id;
    dfloat *ze = mesh->EZ + id;

    for(iint k=0;k<gjNq;++k){
      for(iint j=0;j<gjNq;++j){
	for(iint i=0;i<gjNq;++i){
	  
	  iint n = i + j*gjNq + k*gjNq*gjNq;

	  /* local node coordinates */
	  dfloat rn = mesh->gjr[i]; 
	  dfloat sn = mesh->gjr[j];
	  dfloat tn = mesh->gjr[k];
	  
	  /* Jacobian matrix */
	  dfloat xr = 0.125*( (1-tn)*(1-sn)*(xe[1]-xe[0]) + (1-tn)*(1+sn)*(xe[2]-xe[3]) + (1+tn)*(1-sn)*(xe[5]-xe[4]) + (1+tn)*(1+sn)*(xe[6]-xe[7]) );
	  dfloat xs = 0.125*( (1-tn)*(1-rn)*(xe[3]-xe[0]) + (1-tn)*(1+rn)*(xe[2]-xe[1]) + (1+tn)*(1-rn)*(xe[7]-xe[4]) + (1+tn)*(1+rn)*(xe[6]-xe[5]) );
	  dfloat xt = 0.125*( (1-rn)*(1-sn)*(xe[4]-xe[0]) + (1+rn)*(1-sn)*(xe[5]-xe[1]) + (1+rn)*(1+sn)*(xe[6]-xe[2]) + (1-rn)*(1+sn)*(xe[7]-xe[3]) );
	  
	  dfloat yr = 0.125*( (1-tn)*(1-sn)*(ye[1]-ye[0]) + (1-tn)*(1+sn)*(ye[2]-ye[3]) + (1+tn)*(1-sn)*(ye[5]-ye[4]) + (1+tn)*(1+sn)*(ye[6]-ye[7]) );
	  dfloat ys = 0.125*( (1-tn)*(1-rn)*(ye[3]-ye[0]) + (1-tn)*(1+rn)*(ye[2]-ye[1]) + (1+tn)*(1-rn)*(ye[7]-ye[4]) + (1+tn)*(1+rn)*(ye[6]-ye[5]) );
	  dfloat yt = 0.125*( (1-rn)*(1-sn)*(ye[4]-ye[0]) + (1+rn)*(1-sn)*(ye[5]-ye[1]) + (1+rn)*(1+sn)*(ye[6]-ye[2]) + (1-rn)*(1+sn)*(ye[7]-ye[3]) );
	  
	  dfloat zr = 0.125*( (1-tn)*(1-sn)*(ze[1]-ze[0]) + (1-tn)*(1+sn)*(ze[2]-ze[3]) + (1+tn)*(1-sn)*(ze[5]-ze[4]) + (1+tn)*(1+sn)*(ze[6]-ze[7]) );
	  dfloat zs = 0.125*( (1-tn)*(1-rn)*(ze[3]-ze[0]) + (1-tn)*(1+rn)*(ze[2]-ze[1]) + (1+tn)*(1-rn)*(ze[7]-ze[4]) + (1+tn)*(1+rn)*(ze[6]-ze[5]) );
	  dfloat zt = 0.125*( (1-rn)*(1-sn)*(ze[4]-ze[0]) + (1+rn)*(1-sn)*(ze[5]-ze[1]) + (1+rn)*(1+sn)*(ze[6]-ze[2]) + (1-rn)*(1+sn)*(ze[7]-ze[3]) );
	  
	  /* compute geometric factors for affine coordinate transform*/
	  dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

	  if(J<1e-12) printf("J = %g !!!!!!!!!!!!!\n", J);
	  
	  dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
	  dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
	  dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
	  
	  dfloat JW = J*mesh->gjw[i]*mesh->gjw[j]*mesh->gjw[k];
	  
	  /* store second order geometric factors */
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
	  gjGeo[NgjGeo*gjNp*e + n + gjNp*GWJID] = JW;
	  
	}
      }
    }
  }

  return gjGeo;
}


void ellipticComputeDegreeVector(mesh3D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatter(mesh, ogs, o_deg, o_deg, dfloatString, "add");
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options) {

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  int NblockV = mymax(1,1024/mesh->Np); // works for CUDA
  int NblockS = mymax(1,1024/maxNodes); // works for CUDA
  int NblockG;

  iint gjNq = mesh->gjNq;
  iint gjNp = gjNq*gjNq*gjNq;
  iint gjNq2 = gjNq*gjNq;
  if(gjNq2<=32) 
    NblockG = ( 32/gjNq2 );
  else {
    if(mesh->Nq<=6) {
      NblockG = 256/gjNq2; 
    }
    else 
      NblockG = 1;
  }
  //  NblockG = 512/gNq2;

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;

  iint Nblock = (Ntotal+blockSize-1)/blockSize;

  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;

  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

  solver->mesh = mesh;

  solver->p   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->r   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->z   = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->zP  = (dfloat*) calloc(NallP, sizeof(dfloat));
  solver->Ax  = (dfloat*) calloc(Nall, sizeof(dfloat));
  solver->Ap  = (dfloat*) calloc(Nall, sizeof(dfloat));

  solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  solver->grad = (dfloat*) calloc(7*(Ntotal+Nhalo), sizeof(dfloat));
  
  solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->r);
  solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
  solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat), solver->zP); // CAUTION
  solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
  solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);
  solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);
  solver->o_pAp  = mesh->device.malloc(sizeof(dfloat));

  solver->o_Aw  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_w    = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  solver->o_s    = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
  
  iint Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  
#if 0
  solver->sendBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
  solver->recvBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
#else
  solver->defaultStream = mesh->device.getStream();
  solver->dataStream = mesh->device.createStream();
  mesh->device.setStream(solver->defaultStream);
  
  if(Nbytes>0){
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
    
    solver->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    solver->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
  }else{
    solver->sendBuffer = NULL;
    solver->recvBuffer = NULL;
  }
  mesh->device.setStream(solver->defaultStream);
#endif
  solver->Nblock = Nblock;

  // BP3 specific stuff starts here
  dfloat *gjGeo = ellipticGeometricFactorsHex3D(mesh);

#if 0
  // build a gjD that maps gjNq to gjNq
  dfloat *gjD2 = (dfloat*) calloc(gjNq*gjNq, sizeof(dfloat));
  for(iint n=0;n<gjNq;++n){
    for(iint m=0;m<gjNq;++m){
      for(iint i=0;i<mesh->Nq;++i){
	gjD2[n*gjNq+m] += mesh->gjD[n*mesh->Nq+i]*mesh->gjI[m*mesh->Nq+i];
      }
    }
  }
#endif

  // TW: temporarily resize gjD
  mesh->gjD = (dfloat*) realloc(mesh->gjD, gjNq*gjNq*sizeof(dfloat)); 
  solver->o_gjD = mesh->device.malloc(gjNq*gjNq*sizeof(dfloat), mesh->gjD);
  solver->o_gjD2 = mesh->device.malloc(gjNq*gjNq*sizeof(dfloat), mesh->gjD2);
  solver->o_gjI = mesh->device.malloc(gjNq*mesh->Nq*sizeof(dfloat), mesh->gjI);
  solver->o_gjGeo = mesh->device.malloc(mesh->Nggeo*gjNp*mesh->Nelements*sizeof(dfloat), gjGeo);
  // BP3 specific stuff ends here 

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  //  kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
  //  kernelInfo.addCompilerFlag("-g");
  kernelInfo.addCompilerFlag("-O3");

  // generically used for blocked DEVICE reductions
  kernelInfo.addDefine("p_blockSize", blockSize);

  kernelInfo.addDefine("p_maxNodes", maxNodes);
  kernelInfo.addDefine("p_Nmax", maxNodes);

  kernelInfo.addDefine("p_NblockV", NblockV);
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_NblockG", NblockG);

  kernelInfo.addDefine("p_Lambda2", 0.5f);

  kernelInfo.addDefine("p_gjNq", mesh->gjNq);
  kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
  kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP*mesh->NqP));
  kernelInfo.addDefine("p_Nverts", mesh->Nverts);

  //  occa::setVerboseCompilation(0);

  for(iint r=0;r<size;++r){
    MPI_Barrier(MPI_COMM_WORLD);
    if(r==rank){
      printf("Building kernels for rank %d\n", rank);
      fflush(stdout);
      mesh->haloExtractKernel =
	saferBuildKernelFromSource(mesh->device, 
				   DHOLMES "/okl/meshHaloExtract3D.okl",
				   "meshHaloExtract3D",
				   kernelInfo);
      
      mesh->gatherKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/gather.okl",
				   "gather",
				   kernelInfo);
      
      mesh->scatterKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/scatter.okl",
				   "scatter",
				   kernelInfo);
      
      mesh->gatherScatterKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/gatherScatter.okl",
				   "gatherScatter",
				   kernelInfo);

      
      mesh->getKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/get.okl",
				   "get",
				   kernelInfo);
      
      mesh->putKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/put.okl",
				   "put",
				   kernelInfo);
      
#if 0
      // WARNING
      if(mesh->Nq<12){
	solver->AxKernel =
	  saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticAxHex3D.okl",
				     "ellipticAxHex3D_e3",
				     kernelInfo);
      }
#endif

      // CPU version is e8, GPU version is e6
      solver->partialAxKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticAxHex3D.okl",
				   "ellipticPartialAxHex3D_e6a", 
				   //				   "ellipticPartialAxHex3D_incremental", 
				   kernelInfo);
      
      
      mesh->weightedInnerProduct1Kernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/weightedInnerProduct1.okl",
				   "weightedInnerProduct1",
				   kernelInfo);
      
      mesh->weightedInnerProduct2Kernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/weightedInnerProduct2.okl",
				   "weightedInnerProduct2",
				   kernelInfo);
      
      mesh->innerProductKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/innerProduct.okl",
				   "innerProduct",
				   kernelInfo);
      
      mesh->scaledAddKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/scaledAdd.okl",
				   "scaledAdd",
				   kernelInfo);
      
      mesh->dotMultiplyKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/dotMultiply.okl",
				   "dotMultiply",
				   kernelInfo);
      
      mesh->dotDivideKernel = 
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/dotDivide.okl",
				   "dotDivide",
				   kernelInfo);
      
      solver->gradientKernel = 
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticGradientHex3D.okl",
				   "ellipticGradientHex3D",
				   kernelInfo);
      
      solver->partialGradientKernel = 
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticGradientHex3D.okl",
				   "ellipticPartialGradientHex3D",
				   kernelInfo);
      
      if(mesh->Nq<12){
	solver->ipdgKernel =
	  saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticAxIpdgHex3D.okl",
					 "ellipticAxIpdgHex3D",
				     kernelInfo);
	
	solver->partialIpdgKernel =
	  saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticAxIpdgHex3D.okl",
				     "ellipticPartialAxIpdgHex3D",
				     kernelInfo);
      }
      
      solver->combinedInnerProductKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticCombinedInnerProduct.okl",
				   "ellipticCombinedInnerProduct",
				   kernelInfo);
      
      solver->combinedUpdateKernel =
	saferBuildKernelFromSource(mesh->device, DHOLMES "/okl/ellipticCombinedUpdate.okl",
				   "ellipticCombinedUpdate",
				   kernelInfo);
      usleep(8000);
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);

  occaTimerTic(mesh->device,"GatherScatterSetup");
  
  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  solver->ogs = meshParallelGatherScatterSetup(mesh,
					       mesh->Np*mesh->Nelements,
					       sizeof(dfloat),
					       mesh->gatherLocalIds,
					       mesh->gatherBaseIds, 
					       mesh->gatherHaloFlags);
  occaTimerToc(mesh->device,"GatherScatterSetup");

  occaTimerTic(mesh->device,"DegreeVectorSetup");
  dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *degree = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  solver->o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);
  
  ellipticComputeDegreeVector(mesh, Ntotal, solver->ogs, degree);

  for(iint n=0;n<Ntotal;++n){ // need to weight inner products{
    if(degree[n] == 0) printf("WARNING!!!!\n");
    invDegree[n] = 1./degree[n];
  }
  
  solver->o_invDegree.copyFrom(invDegree);
  occaTimerToc(mesh->device,"DegreeVectorSetup");
  
  //fill geometric factors in halo
  if(mesh->totalHaloPairs){
    iint Nlocal = mesh->Nelements*mesh->Np;
    iint Nhalo = mesh->totalHaloPairs*mesh->Np;

    dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));
    
    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));
    
    meshHaloExchange(mesh,
         mesh->Nvgeo*mesh->Np*sizeof(dfloat),
         mesh->vgeo,
         vgeoSendBuffer,
         mesh->vgeo + Nlocal*mesh->Nvgeo);
    
    mesh->o_vgeo =
      mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
  }

  // build weights for continuous SEM L2 project --->
  dfloat *localMM = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat wJ = mesh->ggeo[e*mesh->Np*mesh->Nggeo + n + GWJID*mesh->Np];
      localMM[n+e*mesh->Np] = wJ;
    }
  }

  occa::memory o_localMM = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  occa::memory o_MM      = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);

  // sum up all contributions at base nodes and scatter back

  ellipticParallelGatherScatter(mesh, solver->ogs, o_localMM, o_MM, dfloatString, "add");

  mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
  mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);

  free(localMM); o_MM.free(); o_localMM.free();

  if(rank==0)
    printf("starting elliptic parallel gather scatter setup\n");

  // set up separate gather scatter infrastructure for halo and non halo nodes
  //  mesh->device.setStream(solver->dataStream);
  ellipticParallelGatherScatterSetup(mesh,
				     mesh->Np*mesh->Nelements,
				     sizeof(dfloat),
				     mesh->gatherLocalIds,
				     mesh->gatherBaseIds, 
				     mesh->gatherHaloFlags,
				     &(solver->halo),
				     &(solver->nonHalo));
  //  mesh->device.setStream(solver->defaultStream);

  
  // count elements that contribute to global C0 gather-scatter
  iint globalCount = 0;
  iint localCount = 0;
  iint *localHaloFlags = (iint*) calloc(mesh->Np*mesh->Nelements, sizeof(int));

  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    localHaloFlags[mesh->gatherLocalIds[n]] += mesh->gatherHaloFlags[n];
  }
  
  for(iint e=0;e<mesh->Nelements;++e){
    iint isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(localHaloFlags[e*mesh->Np+n]>0){
	isHalo = 1;
      }
      if(localHaloFlags[e*mesh->Np+n]<0){
	printf("found halo flag %d\n", localHaloFlags[e*mesh->Np+n]);
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }
  
  //  printf("local = %d, global = %d\n", localCount, globalCount);
  
  solver->globalGatherElementList    = (iint*) calloc(globalCount, sizeof(iint));
  solver->localGatherElementList = (iint*) calloc(localCount, sizeof(iint));
  
  globalCount = 0;
  localCount = 0;
  
  for(iint e=0;e<mesh->Nelements;++e){
    iint isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(localHaloFlags[e*mesh->Np+n]>0){
	isHalo = 1;
      }
    }
    if(isHalo){
      solver->globalGatherElementList[globalCount++] = e;
    }
    else{
      solver->localGatherElementList[localCount++] = e;
    }
  }
  //  printf("local = %d, global = %d\n", localCount, globalCount);
  
  solver->NglobalGatherElements = globalCount;
  solver->NlocalGatherElements = localCount;

  if(globalCount)
    solver->o_globalGatherElementList =
      mesh->device.malloc(globalCount*sizeof(iint), solver->globalGatherElementList);
  
  if(localCount)
    solver->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(iint), solver->localGatherElementList);
  
  free(localHaloFlags);
  
  return solver;
}