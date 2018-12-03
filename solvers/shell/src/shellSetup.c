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

#include "shell.h"
#include "omp.h"
#include <unistd.h>

static void shellSetupRHSASBF(shell_t *shell);
static void shellSetupRHSSEM(shell_t *shell);

static dfloat shellManufacturedForcingFunction(shell_t *shell, dfloat x, dfloat y, dfloat z);
static dfloat shellManufacturedForcingFunctionDD(shell_t *shell, dfloat x, dfloat y, dfloat z);
static dfloat shellManufacturedForcingFunctionNN(shell_t *shell, dfloat x, dfloat y, dfloat z);

shell_t *shellSetup(mesh_t *mesh, dfloat lambda, occa::properties kernelInfo, setupAide options)
{
  shell_t *shell = new shell_t();
  shell->mesh = mesh;
  shell->options = options;

  options.getArgs("MESH DIMENSION", shell->dim);
  options.getArgs("ELEMENT TYPE", shell->elementType);

  mesh->Nfields = 1;

  if(shell->dim==3){
    if(shell->elementType == QUADRILATERALS){
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
    }else if(shell->elementType == TRIANGLES){
      meshOccaSetupTri3D(mesh, options, kernelInfo);
    }else{
      meshOccaSetup3D(mesh, options, kernelInfo);
    }
  }
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  occa::streamTag setupStart = mesh->device.tagStream();

  // Set up radial mesh.
  dfloat R;
  options.getArgs("OUTER RADIUS", R);

  if (options.compareArgs("SHELL SOLVER", "SEM") ||
      (options.compareArgs("SHELL SOLVER", "ASBF") &&
       options.compareArgs("ASBF RADIAL BASIS", "PIECEWISEDISCRETE"))) {
    options.getArgs("RADIAL DISCRETIZATION SIZE", shell->Nradelements);
    if (options.compareArgs("RADIAL DISCRETIZATION", "EQUISPACED")) {
      dfloat h = (R - 1.0)/shell->Nradelements;
      shell->Rbreaks = (dfloat*)calloc(shell->Nradelements + 1, sizeof(dfloat));
      shell->Rbreaks[0] = 1.0;
      for (int i = 1; i < shell->Nradelements; i++)
        shell->Rbreaks[i] = 1.0 + i*h;
      shell->Rbreaks[shell->Nradelements] = R;
    } else {
      printf("ERROR:  Unrecognized value for RADIAL DISCRETIZATION option.\n");
      exit(-1);
    }
  } else {
    shell->Nradelements = 0;
    shell->Rbreaks = NULL;
  }

  // Boundary conditions.  (1 for Dirichlet, 2 for Neumann)
  if (options.compareArgs("INNER BC", "DIRICHLET")) {
    shell->innerBC = 1;
  } else if (options.compareArgs("INNER BC", "NEUMANN")) {
    shell->innerBC = 2;
  } else {
    printf("ERROR:  Invalid value for INNER BC option.\n");
    exit(-1);
  }

  if (options.compareArgs("OUTER BC", "DIRICHLET")) {
    shell->outerBC = 1;
  } else if (options.compareArgs("OUTER BC", "NEUMANN")) {
    shell->outerBC = 2;
  } else {
    printf("ERROR:  Invalid value for OUTER BC option.\n");
    exit(-1);
  }

  // Set up the solver.
  occa::streamTag solveSetupStart = mesh->device.tagStream();
  shellSolveSetup(shell, lambda, kernelInfo);
  occa::streamTag solveSetupEnd = mesh->device.tagStream();
  shell->times.setup.solveSetup = mesh->device.timeBetween(solveSetupStart, solveSetupEnd);

  // Set up the right-hand side.
  occa::streamTag rhsSetupStart = mesh->device.tagStream();
  if (options.compareArgs("SHELL SOLVER", "ASBF")) {
    shellSetupRHSASBF(shell);
  } else if (options.compareArgs("SHELL SOLVER", "SEM")) {
    shellSetupRHSSEM(shell);
  } else {
    printf("ERROR:  Invalid value \"%s\" for SHELL SOLVER.\n",
           options.getArgs("SHELL SOLVER").c_str());
    exit(-1);
  }
  occa::streamTag rhsSetupEnd = mesh->device.tagStream();
  shell->times.setup.rhsSetup = mesh->device.timeBetween(rhsSetupStart, rhsSetupEnd);

  occa::streamTag setupEnd = mesh->device.tagStream();
  shell->times.setup.total = mesh->device.timeBetween(setupStart, setupEnd);

  return shell;
}

/******************************************************************************/

static void shellSetupRHSASBF(shell_t *shell)
{
  mesh_t *mesh = shell->mesh;

  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat xbase = mesh->x[e*mesh->Np+n];
      dfloat ybase = mesh->y[e*mesh->Np+n];
      dfloat zbase = mesh->z[e*mesh->Np+n];

      dfloat J;

      if(shell->elementType==QUADRILATERALS)
        J = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JID) + n];
      else
        J = mesh->vgeo[e*mesh->Nvgeo + JID];

      for(int g=0;g<shell->Nquad;++g){

        dfloat Rg = shell->Rquad[g];

        // stretch coordinates
        dfloat xg = Rg*xbase;
        dfloat yg = Rg*ybase;
        dfloat zg = Rg*zbase;

        // evaluate rhs at shell quadrature for each surface node
        shell->f[g] = shellManufacturedForcingFunction(shell, xg, yg, zg);
      }

      // integrate f against shell modes
      for(int m=0;m<shell->Nmodes;++m){
        dfloat fhatm = 0;
        for(int i=0;i<shell->Nquad;++i){
          fhatm += shell->Bquad[m + i*shell->Nmodes]*shell->Wquad[i]*shell->f[i];
        }

        // scale by surface weight
        shell->r3D[e*mesh->Np + n + m*shell->Ntotal] = J*fhatm;
      }
    }
  }

  return;
}

static void shellSetupRHSSEM(shell_t *shell)
{
  elliptic_t *elliptic = shell->elliptic;
  mesh_t *mesh         = elliptic->mesh;

  for (int e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) {
      dfloat J = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JID) + n];
      
      dlong ind = e*mesh->Np + n;
      dfloat xn = mesh->x[ind];
      dfloat yn = mesh->y[ind];
      dfloat zn = mesh->z[ind];

      elliptic->r[ind] = J*shellManufacturedForcingFunction(shell, xn, yn, zn);
    }
  }
}

/*****************************************************************************/

// NB:  This must match shellManufacturedSolution() in shellErrorHex3D.c
static dfloat shellManufacturedForcingFunction(shell_t *shell, dfloat x, dfloat y, dfloat z)
{
  if ((shell->innerBC == 1) && (shell->outerBC == 1)) {
    return shellManufacturedForcingFunctionDD(shell, x, y, z);
  } else if ((shell->innerBC == 2) && (shell->outerBC == 2)) {
    return shellManufacturedForcingFunctionNN(shell, x, y, z);
  } else {
    printf("ERROR:  No manufactured solution for type-(%d, %d) BCs.\n",
           shell->innerBC, shell->outerBC);
    exit(-1);
  }
}

// Appropriate BCs:  Dirichlet-Dirichlet
static dfloat shellManufacturedForcingFunctionDD(shell_t *shell, dfloat x, dfloat y, dfloat z)
{
  dfloat r, theta, phi;
  dfloat f;

  r     = sqrt(x*x + y*y + z*z);
  theta = atan2(y, x);
  phi   = acos(z/r);

#if 0
  // NB:  This solution assumes shell->R == 1.5.
  dfloat k1 = 6.283185307179586;
  dfloat k2 = 18.849555921538759;
  dfloat k3 = 25.132741228718345;
  f = (k1 + shell->lambda/k1)*sin(k1*r)/r
          + (k2 + shell->lambda/k2)*sin(k2*r)/r
          + (k3 + shell->lambda/k3)*sin(k3*r)/r;
#elif 1
  dfloat q, d2qdx2, d2qdy2, d2qdz2;
  dfloat r2mR2      = pow(r, 2.0) - pow(shell->R, 2.0);
  dfloat r2m1       = pow(r, 2.0) - 1.0;
  dfloat twor2mR2m1 = 2.0*pow(r, 2.0) - pow(shell->R, 2.0) - 1.0;
  
  q = sin(x)*cos(y)*exp(z)*r2mR2*r2m1;
  d2qdx2 = cos(y)*exp(z)*((2.0*sin(x) + 4.0*x*cos(x))*twor2mR2m1 + 8.0*x*x*sin(x) - sin(x)*r2mR2*r2m1);
  d2qdy2 = sin(x)*exp(z)*((2.0*cos(y) - 4.0*y*sin(y))*twor2mR2m1 + 8.0*y*y*cos(y) - cos(y)*r2mR2*r2m1);
  d2qdz2 = sin(x)*cos(y)*exp(z)*((2.0 + 4.0*z)*twor2mR2m1 + 8.0*z*z + r2mR2*r2m1);
  f = -d2qdx2 - d2qdy2 - d2qdz2 + shell->lambda*q;
#endif

  return f;
}

// Appropriate BCs:  Neumann-Neumann
static dfloat shellManufacturedForcingFunctionNN(shell_t *shell, dfloat x, dfloat y, dfloat z)
{
  dfloat r, theta, phi;
  dfloat f;

  r     = sqrt(x*x + y*y + z*z);
  theta = atan2(y, x);
  phi   = acos(z/r);

  // NB:  This solution may only work for shell->R not too much bigger than 1.
  dfloat q, dqdr, d2qdr2, dqdphi, d2qdphi2, d2qdtheta2;
  dfloat lapq, lapqr, lapqtheta, lapqphi;
  dfloat p, dpdr, d2pdr2;
  dfloat s, dsdphi, d2sdphi2;

#if 0
  p = r*(pow(r, 2.0)/3.0 - ((1.0 + shell->R)/2.0)*r + shell->R);
  dpdr = (1.0 - r)*(shell->R - r);
  d2pdr2 = 2.0*r - (1.0 + shell->R);
#else
  p = cos(M_PI*(r - 1.0)/(shell->R - 1.0));
  dpdr = -M_PI*sin(M_PI*(r - 1.0)/(shell->R - 1.0))/(shell->R - 1.0);
  d2pdr2 = -M_PI*M_PI*cos(M_PI*(r - 1.0)/(shell->R - 1.0))/pow(shell->R - 1.0, 2.0);
#endif

  s = pow(phi, 3.0)*pow(phi - M_PI, 3.0);
  dsdphi = 3.0*pow(phi, 2.0)*pow(phi - M_PI, 2.0)*(2.0*phi - M_PI);
  d2sdphi2 = 6.0*phi*(phi - M_PI)*(5.0*pow(phi, 2.0) - 5.0*phi*M_PI + M_PI*M_PI);

  q = sin(theta)*s*p;
  dqdr = sin(theta)*s*dpdr;
  d2qdr2 = sin(theta)*s*d2pdr2;
  d2qdtheta2 = -sin(theta)*s*p;
  dqdphi = sin(theta)*dsdphi*p;
  d2qdphi2 = sin(theta)*d2sdphi2*p;

  if ((fabs(phi) < 1.0e-13) || (fabs(phi - M_PI) < 1.0e-13)) {
    // Near poles, the Laplacian is approximately zero.
    lapq = 0;
  } else {
    lapqr = d2qdr2 + (2.0/r)*dqdr;
    lapqphi = (1.0/pow(r, 2.0))*d2qdphi2 + (cos(phi)/(pow(r, 2.0)*sin(phi)))*dqdphi;
    lapqtheta = d2qdtheta2/(pow(r, 2.0)*pow(sin(phi), 2.0));
    lapq = lapqr + lapqtheta + lapqphi;
  }

  f = -lapq + shell->lambda*q;

  return f;
}
