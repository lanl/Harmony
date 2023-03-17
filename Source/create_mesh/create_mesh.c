/**********************************************************************************************
Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and 
to permit others to do so.
**********************************************************************************************/
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define INDEX_FOR_3D(x,xl,xh,y,yl,yh,z,zl) (((xh)-(xl)+1)*((y)-(yl))+(x)-(xl) + ((xh)-(xl)+1)*((yh)-(yl)+1)*((z)-(zl)))

typedef double realtyp;

typedef int64_t INT;

#include "parser.h"


void write_ensight(char *dir,char *name,int num_nodes, int num_elements,double *x, double *y, double *z,
                    int *connect)
{

  int i,j,k,l,id;

  FILE *out[20];

  double xmin,xmax;

  char file[100];

  sprintf(file,"%s%s.geo",dir,name);
  
  out[0]=fopen(file,"w");


  fprintf(out[0],"EnSight Model Geometry File\n");
  fprintf(out[0],"EnSight 10.1.6\n");
  fprintf(out[0],"node id given\n");
  fprintf(out[0],"element id given\n");
  fprintf(out[0],"extents\n");

  xmin=1.e5;
  xmax=-1.e5;
  for (id=0;id<num_nodes;id++) {
    if (x[id]<xmin) xmin=x[id];
    if (x[id]>xmax) xmax=x[id];
  }
    
  fprintf(out[0],"%12.5e%12.5e\n",xmin,xmax);
  xmin=1.e5;
  xmax=-1.e5;
  for (id=0;id<num_nodes;id++) {
    if (y[id]<xmin) xmin=y[id];
    if (y[id]>xmax) xmax=y[id];
  }
    
  fprintf(out[0],"%12.5e%12.5e\n",xmin,xmax);
  xmin=1.e5;
  xmax=-1.e5;
  for (id=0;id<num_nodes;id++) {
    if (z[id]<xmin) xmin=z[id];
    if (z[id]>xmax) xmax=z[id];
  }
    
  fprintf(out[0],"%12.5e%12.5e\n",xmin,xmax);
  
  
  fprintf(out[0],"part\n");
  fprintf(out[0],"%10d\n",1);
  fprintf(out[0],"Block ID 100000 - HEX\n");
  fprintf(out[0],"coordinates\n");
  fprintf(out[0],"%10lu\n",num_nodes);

  for (id=1;id<=num_nodes;id++)
    fprintf(out[0],"%10lu\n",id);
  
  for (id=0;id<num_nodes;id++)
    fprintf(out[0],"%12.5e\n",x[id]);

  for (id=0;id<num_nodes;id++)
    fprintf(out[0],"%12.5e\n",y[id]);
  
  for (id=0;id<num_nodes;id++)
    fprintf(out[0],"%12.5e\n",z[id]);  


  fprintf(out[0],"hexa8\n");
  fprintf(out[0],"%10lu\n",num_elements);

  for (id=1;id<=num_elements;id++)
    fprintf(out[0],"%10lu\n",id);  

  k=0;
  for (id=0;id<num_elements;id++){

    fprintf(out[0],"%10lu",connect[k]);
    fprintf(out[0],"%10lu",connect[k+1]);
    fprintf(out[0],"%10lu",connect[k+2]);
    fprintf(out[0],"%10lu",connect[k+3]);

    fprintf(out[0],"%10lu",connect[k+4]);
    fprintf(out[0],"%10lu",connect[k+5]);
    fprintf(out[0],"%10lu",connect[k+6]);
    fprintf(out[0],"%10lu",connect[k+7]);
    fprintf(out[0],"\n");

    k=k+8;
   
  }
  fclose(out[0]);
  
  sprintf(file,"%s%s.case",dir,name);
  
  out[0]=fopen(file,"w");

  fprintf(out[0],"# Wed Aug 29 23:50:44 2018\n");
  fprintf(out[0],"# EnSight Gold Model\n");
  fprintf(out[0],"# Produced with EnSight Ouput API, version 1.0.4\n");
  fprintf(out[0],"\n");
  fprintf(out[0],"FORMAT\n");
  fprintf(out[0],"type:                   ensight gold\n");
  fprintf(out[0],"\n");
  fprintf(out[0],"GEOMETRY\n");
  fprintf(out[0],"model:                  %s.geo\n",name);
  //  fprintf(out[0],"VARIABLE\n");
  //  fprintf(out[0],"scalar per element: flg %s.flg\n",name);
  fclose(out[0]);
  
}


void partition(int num_nodes, int num_elements,double *x, double *y, double *z,int *connect)
{

  int i,j,k,l,ibkg,ii,jj,kk;
  int i1,i2,j1,j2,k1,k2;
  int is1,is2,js1,js2,ks1,ks2;  

  int dx,dy,dz;

  int num_cells;

  int num_points;

  int num;

  int *cell;

  int *point;

  int *flag;

  size_t size;

  int rank,id,m,gid;

  double *x_local,*y_local,*z_local;

  int *connect_local;

  int *global_to_local;

  char dir[100],name[100];

  FILE *out;

  global_to_local=(int *) malloc((size_t) (num_nodes*sizeof(int)));
  
  dx=gd.LIMX/gd.NPEX;
  dy=gd.LIMY/gd.NPEY;  
  dz=gd.LIMZ/gd.NPEZ;

  rank=0;

  for (k=0;k<gd.NPEZ;k++)
    for (j=0;j<gd.NPEY;j++)
      for (i=0;i<gd.NPEX;i++) {

	//indices for owned cells
	is1=i1=i*dx;
	is2=i2=(i+1)*dx;
	
	js1=j1=j*dy;
	js2=j2=(j+1)*dy;

	ks1=k1=k*dz;
	ks2=k2=(k+1)*dz;

	//indices for owned and shared cells
        if (i1!=0) is1=i1-1;
	if (i2!=gd.LIMX) is2=i2+1;
	
        if (j1!=0) js1=j1-1;
	if (j2!=gd.LIMY) js2=j2+1;
	
        if (k1!=0) ks1=k1-1;
	if (k2!=gd.LIMZ) ks2=k2+1;

        	
	num_cells=(is2-is1)*(js2-js1)*(ks2-ks1);
	num_points=(is2-is1+1)*(js2-js1+1)*(ks2-ks1+1);
	
        flag=(int *) malloc((size_t) (num_cells*sizeof(int)));
        connect_local=(int *) malloc((size_t) (8*num_cells*sizeof(int)));
	
	x_local=(double *) malloc((size_t) (num_points*sizeof(double)));
	y_local=(double *) malloc((size_t) (num_points*sizeof(double)));
	z_local=(double *) malloc((size_t) (num_points*sizeof(double)));

	//create ensight directory for this rank
	sprintf(name,"ls ensight_%05d",rank);
	l=system(name);

	if (l!=0) {
	  sprintf(name,"mkdir ensight_%05d",rank);
	  l=system(name);
	}
	
	//labels local points by global id
	for(kk=0;kk<num_nodes;kk++)
	  global_to_local[kk]=-1;

	sprintf(dir,"ensight_%05d/engold.gid",rank);
	out=fopen(dir,"w"); //write file with list of global ids 
	
	id=0;
	for (kk=ks1;kk<=ks2;kk++)
	  for (jj=js1;jj<=js2;jj++)
	    for (ii=is1;ii<=is2;ii++) {      
	      l=INDEX_FOR_3D(ii,0,gd.LIMX,jj,0,gd.LIMY,kk,0);
	      global_to_local[l]=id;
	      x_local[id]=x[l];
	      y_local[id]=y[l];
	      z_local[id]=z[l];
	      fprintf(out,"%d\n",l);
	      id++;
	    }
	fclose(out);
	printf("rank= %d NPoints= %d id= %d\n",rank,num_points,id);	
	
	//labels local cells by global id
	sprintf(dir,"ensight_%05d/engold.cid",rank);
	out=fopen(dir,"w"); //write file with list of global ids 
	
	id=0;
	for (kk=ks1;kk<ks2;kk++)
	  for (jj=js1;jj<js2;jj++)
	    for (ii=is1;ii<is2;ii++) {
	      
	      gid=INDEX_FOR_3D(ii,0,gd.LIMX-1,jj,0,gd.LIMY-1,kk,0);
	      fprintf(out,"%d\n",gid);
	      
	      flag[id]=0;
	      if (kk>=k1 && kk<k2 && jj>=j1 && jj<j2 && ii>=i1 && ii<i2) flag[id]=1;
	      
	      for (l=0;l<8;l++) {
		m=connect[8*gid+l]-1;
		connect_local[8*id+l]=global_to_local[m]+1;
	      }
	      id++;
	    }
	fclose(out);
	printf("rank= %d NCells= %d id= %d\n",rank,num_cells,id);
	
        //write ensight geometry file for this rank
	sprintf(dir,"ensight_%05d/",rank);
	sprintf(name,"engold");	
	write_ensight(dir,name,num_points,num_cells,x_local,y_local,z_local,connect_local);

        //write cell ownership flag
	sprintf(dir,"ensight_%05d/engold.flg",rank);
	out=fopen(dir,"w"); //write file with list of global ids
	//fprintf(out,"Per_elem scalar values\n");
	//fprintf(out,"part\n");
	//fprintf(out,"%10d\n",1);
	//fprintf(out,"hexa8\n");
	
	id=0;
	for (kk=ks1;kk<ks2;kk++)
	  for (jj=js1;jj<js2;jj++)
	    for (ii=is1;ii<is2;ii++) {
	      fprintf(out,"%d\n",flag[id]);
	      id++;
	    }
	fclose(out);
	      
	
	free(flag);
	free(connect_local);	
	free(x_local);
	free(y_local);
	free(z_local);

	rank++;
      }//end loop over ranks

  free(global_to_local);

}

int set_mesh()
{

  int i,j,k,l,ibkg,ii,jj,kk;
  int i1,i2,j1,j2,k1,k2,JT;

  double dx,dy,dz;
  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  double rho,E,mass;

  double theta1, theta2, theta, dtheta;
  double phi1, phi2, phi, dphi;
  double avex,avey,avez;
 
  double *th1,*r1,*yp1,rad,Pi;
  int num;

  double *x,*y,*z;

  int cnt,debug;

  int *connect;

  int num_elements;
  int num_nodes;
  size_t size;

  char dir[100];
  char name[100];

  sprintf(name,"mesh");
  sprintf(dir,"./");  

  num_nodes=(gd.LIMX+1)*(gd.LIMY+1)*(gd.LIMZ+1);

  x=(double *) malloc((size_t) (num_nodes*sizeof(double)));
  y=(double *) malloc((size_t) (num_nodes*sizeof(double)));
  z=(double *) malloc((size_t) (num_nodes*sizeof(double)));
  assert(x != NULL && y != NULL && z != NULL);

  Pi=acos(0.)/90.0;

  for (l=0;l<fd.NF;l++) {

    printf("l= %d\n",l);
    
    if (fd.fill[l].type==10) { //cylindrical geometry

      i1=fd.fill[l].i1;
      i2=fd.fill[l].i2;

      j1=fd.fill[l].j1;
      j2=fd.fill[l].j2;

      k1=fd.fill[l].k1;
      k2=fd.fill[l].k2;

      x1=fd.fill[l].theta1*Pi;
      x2=fd.fill[l].theta2*Pi;
      dx=(x2-x1)/(i2-i1);

      y1=fd.fill[l].rad1;
      y2=fd.fill[l].rad2;
      dy=(y2-y1)/(j2-j1);

      z1=fd.fill[l].z1;
      z2=fd.fill[l].z2;
      dz=(z2-z1)/(k2-k1);

      printf("dx= %f dy= %f dz= %f\n",dx,dy,dz);

      for (k=k1;k<=k2;k++)
	for (j=j1;j<=j2;j++) 
          for (i=i1;i<=i2;i++) {
            
            num=INDEX_FOR_3D(i,0,gd.LIMX,j,0,gd.LIMY,k,0);  

	    theta1=x1+dx*(i-i1);
	    rad=y1+dy*(j-j1);
	    x[num]=rad*cos(theta1);
	    y[num]=rad*sin(theta1);
	    z[num]=z1+dz*(k-k1);

	    //printf("%d %d %d\n",k,j,i);

	}

    } //end if type==0
    else if (fd.fill[l].type==1) {

      i1=fd.fill[l].i1;
      i2=fd.fill[l].i2;

      j1=fd.fill[l].j1;
      j2=fd.fill[l].j2;

      k1=fd.fill[l].k1;
      k2=fd.fill[l].k2;

      x1=fd.fill[l].x1;
      x2=fd.fill[l].x2;
      dx=(x2-x1)/(i2-i1);

      y1=fd.fill[l].y1;
      y2=fd.fill[l].y2;
      dy=(y2-y1)/(j2-j1);

      z1=fd.fill[l].z1;
      z2=fd.fill[l].z2;
      dz=(z2-z1)/(k2-k1);

      //printf("dx= %f dy= %f dz= %f\n",dx,dy,dz);

      for (k=k1;k<=k2;k++)
	for (j=j1;j<=j2;j++) 
          for (i=i1;i<=i2;i++) {
            
            num=INDEX_FOR_3D(i,0,gd.LIMX,j,0,gd.LIMY,k,0);  

	    x[num]=x1+dx*(i-i1);
	    y[num]=y1+dy*(j-j1);
	    z[num]=z1+dz*(k-k1);

	}

    } //end if type==1

  } //end NF loop

  num_elements = (gd.LIMX)*(gd.LIMY)*(gd.LIMZ);
  size         = (size_t)8 * num_elements * sizeof(int);
  assert(size > 0);
  connect = malloc(size);
  assert(connect != NULL);

  cnt = 0;
  for (k=1;k<=gd.LIMZ;k++)
    for (j=1;j<=gd.LIMY;j++) 
      for (i=1;i<=gd.LIMX;i++) {

        connect[cnt++] = 1+INDEX_FOR_3D(i-1,0,gd.LIMX,j-1,0,gd.LIMY,k-1,0);  
        connect[cnt++] = 1+INDEX_FOR_3D(i,0,gd.LIMX,j-1,0,gd.LIMY,k-1,0);
        connect[cnt++] = 1+INDEX_FOR_3D(i,0,gd.LIMX,j,0,gd.LIMY,k-1,0);
        connect[cnt++] = 1+INDEX_FOR_3D(i-1,0,gd.LIMX,j,0,gd.LIMY,k-1,0);

        connect[cnt++] = 1+INDEX_FOR_3D(i-1,0,gd.LIMX,j-1,0,gd.LIMY,k,0);  
        connect[cnt++] = 1+INDEX_FOR_3D(i,0,gd.LIMX,j-1,0,gd.LIMY,k,0);
        connect[cnt++] = 1+INDEX_FOR_3D(i,0,gd.LIMX,j,0,gd.LIMY,k,0);
        connect[cnt++] = 1+INDEX_FOR_3D(i-1,0,gd.LIMX,j,0,gd.LIMY,k,0);

      }

  /*
   *    Write Out Mesh
   */

  write_ensight(dir,name,num_nodes,num_elements,x,y,z,connect);
  partition(num_nodes,num_elements,x,y,z,connect);
  

  free(x);
  free(y);
  free(z);
  free(connect);  

  return 0;

}



int main(int argc, char** argv) 
{

  char file[100];



  //read input deck
  parser(argv[1]);

  //write exodus file
  sprintf(file,"my_mesh");
  printf("my_mesh\n");
  set_mesh();

  return 0;

}
