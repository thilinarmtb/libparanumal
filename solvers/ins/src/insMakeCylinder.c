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

void insMakeCylinder(ins_t *ins, setupAide options){
  //
  dfloat xc, yc, R;

  mesh_t *mesh = ins->mesh; 

  int curvedBcId = 1; // change it later

  // Currently Only Quad and Simple Cylinder Setup
  if(ins->elementType==QUADRILATERALS){

    options.getArgs("CYLINDER CENTER X", xc); 
    options.getArgs("CYLINDER CENTER Y", yc); 
    options.getArgs("CYLINDER RADIOUS", R ); 


    dfloat *fdx = (dfloat *) calloc(mesh->Nfp, sizeof(dfloat));
    dfloat *fdy = (dfloat *) calloc(mesh->Nfp, sizeof(dfloat));
     
    printf("Cylinder Radious:%.2e Center X: %.2e Center Y: %.2e\n", R, xc, yc);
    // Find number of Elements needed to be deformed
    // In proper mesh there is only one boundary!!!!
    int cylNfields = 2; // currently element and face ids
    int cylNelements = 0; 
    for(dlong e=0; e<mesh->Nelements; e++){
      for(int f=0; f<mesh->Nfaces; f++){
        if(mesh->EToB[e*mesh->Nfaces + f] == curvedBcId)
          cylNelements++;
      }
    }

    printf("Number of elements having curved boundary: %d\n", cylNelements);

    // Allocate memory for elements and face ids
    int cnt = 0; 
    int *cylIds = (int *) calloc(cylNfields*cylNelements, sizeof(int));
    for(dlong e=0; e<mesh->Nelements; e++){
      for(int f=0; f<mesh->Nfaces; f++){
         if(mesh->EToB[e*mesh->Nfaces + f] == curvedBcId){
          cylIds[cnt++] = e;
          cylIds[cnt++] = f;
         }
       }
     }

     printf("Number of elements having curved boundary: %d\n", cylNelements);

     // Deform face and then blend to volume nodes
     for(int elm=0; elm<cylNelements; elm++){
     // for(int elm=0; elm<1; elm++){
      const int e = cylIds[elm*cylNfields + 0]; 
      const int f = cylIds[elm*cylNfields + 1];

      // printf("Elm :%d  and Face: %d \n", e,f);

      dfloat  x1 = 0.0, x2 = 0.0; 
      dfloat  y1 = 0.0, y2 = 0.0; 
      if(f==0){
        x1 = mesh->EX[e*mesh->Nverts + 0];
        x2 = mesh->EX[e*mesh->Nverts + 1];
        y1 = mesh->EY[e*mesh->Nverts + 0];
        y2 = mesh->EY[e*mesh->Nverts + 1];
      }
      else if(f==1){
        x1 = mesh->EX[e*mesh->Nverts + 1];
        x2 = mesh->EX[e*mesh->Nverts + 2];
        y1 = mesh->EY[e*mesh->Nverts + 1];
        y2 = mesh->EY[e*mesh->Nverts + 2];
      } 
      else if(f==2){
        x1 = mesh->EX[e*mesh->Nverts + 3];
        x2 = mesh->EX[e*mesh->Nverts + 2];
        y1 = mesh->EY[e*mesh->Nverts + 3];
        y2 = mesh->EY[e*mesh->Nverts + 2];
      } 
      else if(f==3){
        x1 = mesh->EX[e*mesh->Nverts + 0];
        x2 = mesh->EX[e*mesh->Nverts + 3];
        y1 = mesh->EY[e*mesh->Nverts + 0];
        y2 = mesh->EY[e*mesh->Nverts + 3];
      }
      
      dfloat theta1 = atan2(y1-yc, x1-xc); 
      dfloat theta2 = atan2(y2-yc, x2-xc);

      dfloat TOL = 1e-10;

      // if ((theta2 > TOL) && (theta1 < TOL)){
      if ((theta2 > TOL) && (theta1 < TOL)){
        printf("there is an quadrand problem on element %d with index %d\n", e, elm);
        printf("theta1 %.4e and theta2 %.4e\n", theta1, theta2);
        printf("theta1 %.4e and theta2 %.4e\n", theta1*180/M_PI, theta2*180/M_PI);
        if(fabs(theta1)>TOL)
         theta1 = theta1 + 2*M_PI;
      }
      if ((theta1 > TOL) && (theta2 < TOL)){
        printf("there is another quadrand problem on element %d with index %d\n", e, elm);
         printf("theta1 %.4e and theta2 %.4e\n", theta1, theta2);
        printf("theta1 %.4e and theta2 %.4e\n", theta1*180/M_PI, theta2*180/M_PI);
        if(fabs(theta2)>TOL)
          theta2 = theta2 + 2*M_PI;

      }

      // printf("Theta1 :%.2e  and Theta: %.2e \n", theta1,theta2);


      for(int n=0; n<mesh->Nfp; n++){
        // compute the radial location of each node on GLL distribution
        dfloat theta = 0.5*theta1*(1.0-mesh->gllz[n]) + 0.5*theta2*(1.0+mesh->gllz[n]);
        int vnode    = mesh->faceNodes[f*mesh->Nfp + n]; 
        fdx[n] = xc + R*cos(theta)- mesh->x[e*mesh->Np + vnode];  
        fdy[n] = yc + R*sin(theta)- mesh->y[e*mesh->Np + vnode];  
      }

      dfloat dx = 0.0, dy = 0.0; 
      // Blend face deformation to volume nodes using linear interpolation
      for(int i=0; i<mesh->Nfp; i++){
        for(int j=0; j<mesh->Nfp; j++){
          const int id = j + i*mesh->Nfp; // r runs faster
          if(f==0){
            dx = 0.5*(1.0-mesh->s[id])*fdx[j];
            dy = 0.5*(1.0-mesh->s[id])*fdy[j];
          }
          else if(f==1){
            dx = 0.5*(1.0+mesh->r[id])*fdx[i];
            dy = 0.5*(1.0+mesh->r[id])*fdy[i];
          } 
          else if(f==2){
            dx = 0.5*(1.0+mesh->s[id])*fdx[j];
            dy = 0.5*(1.0+mesh->s[id])*fdy[j];
          } 
          else if(f==3){
            dx = 0.5*(1.0-mesh->r[id])*fdx[i];
            dy = 0.5*(1.0-mesh->r[id])*fdy[i];
          }

          mesh->x[e*mesh->Np + id] += dx;
          mesh->y[e*mesh->Np + id] += dy;
        }
      }
    }


  
    // Update volume geometric factors
    for(int elm=0; elm<cylNelements; elm++){
     // for(int elm=0; elm<1; elm++){
      const int e = cylIds[elm*cylNfields + 0]; 
      const int f = cylIds[elm*cylNfields + 1];
      // for each node on the element
      for(int n=0; n<mesh->Np; n++){
        dfloat xr = 0.0, xs = 0.0;  
        dfloat yr = 0.0, ys = 0.0;  
        for(int m=0; m<mesh->Np; m++){
          xr +=mesh->Dr[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          xs +=mesh->Ds[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          yr +=mesh->Dr[n*mesh->Np +m]*mesh->y[e*mesh->Np + m]; 
          ys +=mesh->Ds[n*mesh->Np +m]*mesh->y[e*mesh->Np + m];  
        }

        dfloat J = xr*ys - xs*yr;
        if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        int i = n%mesh->Nq;
        int j = n/mesh->Nq;
        dfloat JW = J*mesh->gllw[i]*mesh->gllw[j];

        // // Sanity Check
        // dfloat rxe  = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] - rx)  ;
        // dfloat rye  = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] - ry)  ;
        // dfloat sxe  = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] - sx)  ;
        // dfloat sye  = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] - sy)  ;
        // dfloat Je   = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID] - J ) ;
        // dfloat JWe  = fabs(mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JWID] -JW) ;
        // dfloat TOL = 1e-8; 
        // if(rxe>TOL)
        //   printf("Difference of rx is %.4e \n", rxe);

        // if(rxe>TOL)
        //   printf("Difference of sx is %.4e \n", sxe);

        // if(rxe>TOL)
        //   printf("Difference of ry is %.4e \n", rye);

        // if(rxe>TOL)
        //   printf("Difference of sy is %.4e \n", sye);

        // if(rxe>TOL)
        //   printf("Difference of J is %.4e \n", Je);

        // if(rxe>TOL)
        //   printf("Difference of JW is %.4e \n", JWe);



        /* store geometric factors */
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] = rx;
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] = ry;
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] = sx;
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] = sy;
        // printf("Jacobian before: %.4e and ",mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]);
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]  = J;
        // printf("after %.4e\n",mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]);
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JWID] = JW;
        mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*IJWID] = 1./JW;

        /* store second order geometric factors */
        mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G00ID] = JW*(rx*rx + ry*ry);
        mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G01ID] = JW*(rx*sx + ry*sy);
        mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G11ID] = JW*(sx*sx + sy*sy);
        mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*GWJID] = JW;
      }
    }



 //  Cubature volume data for deformed elements
     dfloat *cxr = (dfloat *) calloc(mesh->Np, sizeof(dfloat));
     dfloat *cxs = (dfloat *) calloc(mesh->Np, sizeof(dfloat));
     dfloat *cyr = (dfloat *) calloc(mesh->Np, sizeof(dfloat));
     dfloat *cys = (dfloat *) calloc(mesh->Np, sizeof(dfloat));

     // hold interpolated variables
     
     dfloat *ixr = (dfloat *) calloc(mesh->cubNq*mesh->Nfp, sizeof(dfloat));
     dfloat *ixs = (dfloat *) calloc(mesh->cubNq*mesh->Nfp, sizeof(dfloat));
     dfloat *iyr = (dfloat *) calloc(mesh->cubNq*mesh->Nfp, sizeof(dfloat));
     dfloat *iys = (dfloat *) calloc(mesh->cubNq*mesh->Nfp, sizeof(dfloat));
     
     
    // Update cubature volume geometric factors
    for(int elm=0; elm<cylNelements; elm++){
      // for(int elm=0; elm<1; elm++){
      const int e = cylIds[elm*cylNfields + 0]; 
      const int f = cylIds[elm*cylNfields + 1];

      for(int n=0; n<mesh->Np; n++){
        cxr[n] = 0.0; cxs[n] = 0.0;  
        cyr[n] = 0.0; cys[n] = 0.0;  
        for(int m=0; m<mesh->Np; m++){
          cxr[n] +=mesh->Dr[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          cxs[n] +=mesh->Ds[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          cyr[n] +=mesh->Dr[n*mesh->Np +m]*mesh->y[e*mesh->Np + m]; 
          cys[n] +=mesh->Ds[n*mesh->Np +m]*mesh->y[e*mesh->Np + m];    
        }
      }
   
    // interpolate in i direction
    for(int j=0;j<mesh->Nfp;++j){
      for(int i=0;i<mesh->cubNq;++i){
        dfloat xr = 0.0, xs = 0.0;  
        dfloat yr = 0.0, ys = 0.0;  
        for (int n=0;n<mesh->Nfp;n++) {
          const int id = j*mesh->Nfp + n; 
          xr += mesh->cubInterp[i*mesh->Nfp+n]*cxr[id];
          xs += mesh->cubInterp[i*mesh->Nfp+n]*cxs[id];
          yr += mesh->cubInterp[i*mesh->Nfp+n]*cyr[id];
          ys += mesh->cubInterp[i*mesh->Nfp+n]*cys[id];
        }
        const int id = j*mesh->cubNq + i; 
        ixr[id] = xr;
        ixs[id] = xs;
        iyr[id] = yr;
        iys[id] = ys;
      }
    }


    for(int j=0; j<mesh->cubNq; j++){
      for(int i=0; i<mesh->cubNq; i++){

        dfloat xr = 0.0, xs = 0.0;  
        dfloat yr = 0.0, ys = 0.0;  
        for (int n=0;n<mesh->Nfp;n++) {
          const int id = n*mesh->cubNq + i; 
          xr += mesh->cubInterp[j*mesh->Nfp+n]*ixr[id];
          xs += mesh->cubInterp[j*mesh->Nfp+n]*ixs[id];
          yr += mesh->cubInterp[j*mesh->Nfp+n]*iyr[id];
          ys += mesh->cubInterp[j*mesh->Nfp+n]*iys[id];
        }

        dfloat J = xr*ys - xs*yr;

        if(J<1e-8) { printf("Negative or small Jacobian: %d %g\n", e, J); exit(-1);}
        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        dfloat JW = J*mesh->cubw[i]*mesh->cubw[j];

         // store geometric factors 
        dlong base = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq;
        
        // // Sanity Check
        // dfloat rxe  = fabs(mesh->cubvgeo[base + mesh->cubNp*RXID] - rx)  ;
        // dfloat rye  = fabs(mesh->cubvgeo[base + mesh->cubNp*RYID] - ry)  ;
        // dfloat sxe  = fabs(mesh->cubvgeo[base + mesh->cubNp*SXID] - sx)  ;
        // dfloat sye  = fabs(mesh->cubvgeo[base + mesh->cubNp*SYID] - sy)  ;
        // dfloat Je   = fabs(mesh->cubvgeo[base + mesh->cubNp*JID]  - J ) ;
        // dfloat JWe  = fabs(mesh->cubvgeo[base + mesh->cubNp*JWID] -JW) ;
        // dfloat TOL = 1e-8; 
        // if(rxe>TOL)
        //   printf("Difference of rx is %.4e \n", rxe);

        // if(sxe>TOL)
        //   printf("Difference of sx is %.4e \n", sxe);

        // if(rye>TOL)
        //   printf("Difference of ry is %.4e \n", rye);

        // if(sye>TOL)
        //   printf("Difference of sy is %.4e \n", sye);

        // if(Je>TOL)
        //   printf("Difference of J is %.4e \n", Je);

        // if(JWe>TOL)
        //   printf("Difference of JW is %.4e \n", JWe);
        
        
        mesh->cubvgeo[base + mesh->cubNp*RXID] = rx;
        mesh->cubvgeo[base + mesh->cubNp*RYID] = ry;
        mesh->cubvgeo[base + mesh->cubNp*SXID] = sx;
        mesh->cubvgeo[base + mesh->cubNp*SYID] = sy;
        // printf("Cubature jacobian before: %.4e and ",mesh->cubvgeo[base + mesh->cubNp*JID]);
        mesh->cubvgeo[base + mesh->cubNp*JID]  = J;
        // printf("after %.4e\n",mesh->cubvgeo[base + mesh->cubNp*JID]);
        mesh->cubvgeo[base + mesh->cubNp*JWID] = JW;
        mesh->cubvgeo[base + mesh->cubNp*IJWID] = 1./JW;
     }
    }


  }



#if 1
     // Update cubature volume geometric factors
    for(int elm=0; elm<cylNelements; elm++){
     // for(int elm=0; elm<1; elm++){
      const int e = cylIds[elm*cylNfields + 0]; 
      const int f = cylIds[elm*cylNfields + 1];

      for(int n=0; n<mesh->Np; n++){
        cxr[n] = 0.0; cxs[n] = 0.0;  
        cyr[n] = 0.0; cys[n] = 0.0;  
        for(int m=0; m<mesh->Np; m++){
          cxr[n] +=mesh->Dr[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          cxs[n] +=mesh->Ds[n*mesh->Np +m]*mesh->x[e*mesh->Np + m]; 
          cyr[n] +=mesh->Dr[n*mesh->Np +m]*mesh->y[e*mesh->Np + m]; 
          cys[n] +=mesh->Ds[n*mesh->Np +m]*mesh->y[e*mesh->Np + m]; 
          
        }
      }

      for(int face=0;face<mesh->Nfaces;++face){ // for each face
        if(face==f){
        // Update face geometric factors
        for(int i=0;i<mesh->Nfp;++i){  // for each node on face

          /* volume index of face node */
          int n = mesh->faceNodes[f*mesh->Nfp+i];
          dfloat xr = cxr[n];
          dfloat xs = cxs[n];
          dfloat yr = cyr[n];
          dfloat ys = cys[n];
          //
          dfloat J = xr*ys - xs*yr;
          dfloat nx = 0.0, ny = 0.0; 
          /* face f normal and length */
          if(f==0){
            nx =  yr; ny = -xr;
          }
          else if(f==1){
            nx =  ys; ny = -xs;
          } 
          else if(f==2){
            nx = -yr; ny =  xr;
          } 
          else if(f==3){
            nx = -ys; ny =  xs;
          }

          dfloat  d = norm2(nx,ny);

          dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        //    // Sanity Check
        // dfloat nxe  = fabs(mesh->sgeo[base + NXID] - nx/d)  ;
        // dfloat nye  = fabs(mesh->sgeo[base + NYID] - ny/d)  ;
        // dfloat sje  = fabs(mesh->sgeo[base + SJID] - d);
        // dfloat ije  = fabs(mesh->sgeo[base + IJID] - 1./J)  ;
       
        // dfloat TOL = 1e-8; 
        // if(nxe>TOL)
        //   printf("Difference of nx is %.4e \n", nxe);

        // if(nye>TOL)
        //   printf("Difference of ny is %.4e \n", nye);

        // if(sje>TOL)
        //   printf("Difference of sj is %.4e \n", sje);

        // if(ije>TOL)
        //   printf("Difference of ij is %.4e \n", ije);

        

          
          mesh->sgeo[base+NXID] = nx/d;
          mesh->sgeo[base+NYID] = ny/d;
          mesh->sgeo[base+SJID] = d;
          mesh->sgeo[base+IJID] = 1./J;

          mesh->sgeo[base+WIJID] = 1./(J*mesh->gllw[0]);
          mesh->sgeo[base+WSJID] = (d/2.)*mesh->gllw[i];

          }

           //geometric data for quadrature
          for(int i=0;i<mesh->cubNq;++i){ 

            dfloat xr = 0.0, xs = 0.0; 
            dfloat yr = 0.0, ys = 0.0; 
            // first interpolate
            for (int j=0;j<mesh->Nfp;j++) {
              /* volume index of face node */
              int n = mesh->faceNodes[f*mesh->Nfp+j];
              xr += mesh->cubInterp[i*mesh->Nfp+j]*cxr[n];
              xs += mesh->cubInterp[i*mesh->Nfp+j]*cxs[n];
              yr += mesh->cubInterp[i*mesh->Nfp+j]*cyr[n];
              ys += mesh->cubInterp[i*mesh->Nfp+j]*cys[n];
            }

            dfloat J = xr*ys - xs*yr;
            dfloat nx = 0.0, ny = 0.0; 
            /* face f normal and length */
            if(f==0){
            nx =  yr; ny = -xr;
            }
            else if(f==1){
            nx =  ys; ny = -xs;
            } 
            else if(f==2){
            nx = -yr; ny =  xr;
            } 
            else if(f==3){
            nx = -ys; ny =  xs;
            }

            dfloat  d = norm2(nx,ny);

            /* output index */
            dlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->cubNq*e + mesh->cubNq*f + i);

            /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
            // printf("nx ny before: %.4e  %.4e and ",mesh->cubsgeo[base+NXID], mesh->cubsgeo[base+NYID]);
            mesh->cubsgeo[base+NXID] = nx/d;
            mesh->cubsgeo[base+NYID] = ny/d;
            // printf("after %.4e %.4e\n",mesh->cubsgeo[base+NXID], mesh->cubsgeo[base+NYID]);
            mesh->cubsgeo[base+SJID] = d;
            // mesh->cubsgeo[base+SJID] = d/2.;
            mesh->cubsgeo[base+IJID] = 1./J;

            mesh->cubsgeo[base+WIJID] = 1./(J*mesh->cubw[0]);
            mesh->cubsgeo[base+WSJID] = (d/2.)*mesh->cubw[i];
          }


        }
      }

    }
#endif
  for(dlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      dlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      dlong baseP = mesh->mapP[baseM];
      if(baseP<0) baseP = baseM;
      
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }  




  

  }  
}