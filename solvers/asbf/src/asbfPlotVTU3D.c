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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "asbf.h"

// interpolate data to plot nodes and save to file (one per process
void asbfPlotVTU3D(asbf_t *asbf, char *fileNameBase, int fld){

  mesh_t *mesh = asbf->mesh;
  
  int rank;
  rank = mesh->rank;

  FILE *fp;
  char fileName[BUFSIZ];
  sprintf(fileName, "%s_%04d.vtu", fileNameBase, rank);
  // strcpy(fileName,fileNameBase);

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");

  int Eloc = (mesh->Nq-1)*(mesh->Nq-1)*(asbf->asbfNplot-1);
  fprintf(fp, "    <Piece NumberOfPoints=\""dlongFormat"\" NumberOfCells=\""dlongFormat"\">\n", 
          mesh->Nelements*mesh->Np*asbf->asbfNplot, 
          mesh->Nelements*Eloc);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh->Nelements;++e){

	// interpolate from asbf to gll nodes
	for(int g=0;g<mesh->Nq;++g){

	  // assume Nfields=1
	  meshSEM->q[e*meshSEM->Np+g*mesh->Np+n] = qg;
	}
      }
    }
    
    for(int p=0;p<asbf->asbfNplot;++p){
      for(int n=0;n<mesh->Np;++n){
	dfloat Rp = asbf->asbfRquad[p];
	
	// stretch coordinates
	dfloat xbase = mesh->x[e*mesh->Np+n];
	dfloat ybase = mesh->y[e*mesh->Np+n];
	dfloat zbase = mesh->z[e*mesh->Np+n];
	
	dfloat xp = Rp*xbase;
	dfloat yp = Rp*ybase;
	dfloat zp = Rp*zbase;
	
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", xp, yp, zp);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int p=0;p<asbf->asbfNplot;++p){
      for(int n=0;n<mesh->Np;++n){
	
	dfloat qp = 0;
	for(int i=0;i<asbf->asbfNmodes;++i){
	  qp += asbf->asbfBplot[i + p*asbf->asbfNmodes]*asbf->q3D[(e*mesh->Np+n)+i*asbf->Ntotal];
	}
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", qp);
      }
    }
  }

  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int k=0;k<asbf->asbfNplot-1;++k){
      for(int j=0;j<mesh->Nq-1;++j){
	for(int i=0;i<mesh->Nq-1;++i){
	  int b = e*mesh->Np*asbf->asbfNplot + i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;
	  fprintf(fp, 
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat"\n ",
		  b,
		  b+1,
		  b+mesh->Nq+1,
		  b+mesh->Nq,
		  b + mesh->Nq*mesh->Nq,
		  b+1 + mesh->Nq*mesh->Nq,
		  b+mesh->Nq+1 + mesh->Nq*mesh->Nq,
		  b+mesh->Nq+ mesh->Nq*mesh->Nq);
	}
      }
    }
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<Eloc;++n){
      cnt += 8;
      fprintf(fp, "       ");
      fprintf(fp, dlongFormat"\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<Eloc;++n){
      fprintf(fp, "12\n"); // HEX code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
