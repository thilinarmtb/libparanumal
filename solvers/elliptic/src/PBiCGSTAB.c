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
#define EPSILON 1e-9

// FROM NEKBONE: not appropriate since it assumes zero initial data
int pbicgstab(elliptic_t* elliptic, dfloat lambda, 
        occa::memory &o_r, occa::memory &o_x, 
        const dfloat tol, const int MAXIT){
  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  
  dfloat norm2;

  /*aux variables */
  occa::memory &o_p      = elliptic->o_p;
  occa::memory &o_phat   = elliptic->o_rtmp;
  occa::memory &o_s      = elliptic->o_s;
  occa::memory &o_shat   = elliptic->o_shat;
  occa::memory &o_rtilde = elliptic->o_z;
  occa::memory &o_t      = elliptic->o_t;
  occa::memory &o_v      = elliptic->o_v;

  occa::memory &o_Ax = elliptic->o_Ax;
  occa::memory &o_invDegree = elliptic->o_invDegree;

  // compute A*x
  ellipticOperator(elliptic,lambda,o_x,elliptic->o_Ax,dfloatString);
  
  // subtract r = b - A*x
  ellipticScaledAdd(elliptic,-1.f,o_Ax,1.f,o_r);

  //set rtilde = r
  o_rtilde.copyFrom(o_r);

  norm2=ellipticWeightedNorm2(elliptic,o_invDegree,o_r);
  dfloat TOL=mymax(tol*tol*norm2,tol*tol);
  if(norm2<TOL) return 0;

  // set scalars rho0,omega0,alpha
  dfloat alpha,beta,rho[2],omega;
  
  int iter;
  for(iter=1;iter<=MAXIT;++iter){
    //rho_i=<rtilde,r_{i-1}>
    rho[0]=ellipticWeightedInnerProduct(elliptic,o_invDegree,o_rtilde,o_r);

    if(fabs(rho[0])<EPSILON) return 2;

    if(iter==1) o_p.copyFrom(o_r);
    else {
      //beta=(rho_i/rho_{i-1})*(alpha/omega_{i-1})
      beta=(rho[0]/rho[1])*(alpha/omega);

      //TODO: Combine
      //p_{i-1}=-omega_{i-1}*v_{i-1}+1.0*p_{i-1}
      ellipticScaledAdd(elliptic,-omega,o_v,1.f,o_p);
      //p_{i}=beta*p_{i-1}+1.0*r_{i-1}
      ellipticScaledAdd(elliptic,1.f,o_r,beta,o_p);
    }
    // Precon.
    //phat=M^{-1}p
    ellipticPreconditioner(elliptic,lambda,o_p,o_phat);
    
    //v_{i}=A*phat
    ellipticOperator(elliptic,lambda,o_phat,o_v,dfloatString);

    //alpha=rho_{i}/<rtilde,v_{i}>
    alpha=rho[0]/ellipticWeightedInnerProduct(elliptic,o_invDegree,o_rtilde,o_v);

    //s=r-alpha*v
    o_s.copyFrom(o_r);
    ellipticScaledAdd(elliptic,-alpha,o_v,1.f,o_s);
    norm2=ellipticWeightedNorm2(elliptic,o_invDegree,o_s);

    if(norm2<TOL){
      //x=x+alpha*phat
      ellipticScaledAdd(elliptic,alpha,o_phat,1.f,o_x);
      return 0;
    }
    //Precon.
    //shat=M^{-1}s
    ellipticPreconditioner(elliptic,lambda,o_s,o_shat);

    //t=A*shat
    ellipticOperator(elliptic,lambda,o_shat,o_t,dfloatString);

    //omega=<t,s>/<t,t>
    omega=ellipticWeightedInnerProduct(elliptic,o_invDegree,o_t,o_s)/ellipticWeightedNorm2(elliptic,o_invDegree,o_t);

    //x+=alpha*phat+beta*shat
    ellipticScaledAdd(elliptic,alpha,o_phat,1.f,o_x);
    ellipticScaledAdd(elliptic,omega,o_shat,1.f,o_x);

    rho[1]=rho[0];
    //r=s-omega*t
    o_r.copyFrom(o_s);
    ellipticScaledAdd(elliptic,-omega,o_t,1.f,o_r);
    norm2=ellipticWeightedNorm2(elliptic,o_invDegree,o_r);

    if(norm2<TOL) return 0;
    if(verbose && mesh->rank==0)
      printf("BiCGStab: it %d r norm %12.12le alpha = %le \n",iter,sqrt(norm2),alpha);    

    if(fabs(omega)<EPSILON) return 3;
  }
  return 1;
}
