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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define INDEX_FOR_1D(x,xl) ((x)-(xl))
#define INDEX_FOR_2D(x,xl,xh,y,yl) (((xh)-(xl)+1)*((y)-(yl))+(x)-(xl))

#define cell(a,b) cell[INDEX_FOR_2D(a,-1,LIMX+2,b,-1)]

typedef struct _cell {
  double p,r,ie;
  double u,v,w;
  double vol;
  int pn[8];
  double xc,yc,zc;
  double flg;
} cell_t;


main(int argc, char *argv[])
{

  int i,j,k,l,id,m,n,p;

  int LIMX,LIMY,LIMZ,JT;

  int kkzln,kkzln2,kkpll,kkpll2;

  FILE *in,*out[100];

  cell_t *cell;

  char name[100],c;

  double pres,xa,ya,za,ua,va,wa;

  double x1,x2,x3,x4,y1,y2,y3,y4;
  double z1,z2,z3,z4;

  double sum,sum1;

  double Pi=2.*acos(0.);

  double *px,*py,*pz;
  double *pu,*pv,*pw;
  double *ip,*jp,*kp;
  double *kfx1,*kfx2,*kfx3;
  double *sxu1,*sxu2,*sxu3;

  int *flg, *pindx;

  int jedit,rank;

  double edtimes[250];

  double etime;

  sprintf(name,"ls ensight");
  i=system(name);

  if (i!=0) {
    sprintf(name,"mkdir ensight");
    i=system(name);
    sprintf(name,"mkdir ensight/data");
    i=system(name);
  }

 

  n=atoi(argv[1]); //number of time edits
  m=atoi(argv[2]); //number of processors
  etime=atof(argv[3]); //graphics dump time increment

  printf("number of time edits= %d\n",n);
  printf("number of processors= %d\n",m);
  printf("edit time increment= %e\n",etime);


  for (jedit=0;jedit<=n;jedit++) {

    sprintf(name,"ensight/data/engold.%05d.geo",jedit);
    out[0]=fopen(name,"w");
    fprintf(out[0],"This is the 1st description line of the EnSight Gold geometry example\n");
    fprintf(out[0],"This is the 2nd description line of the EnSight Gold geometry example\n");
    fprintf(out[0],"node id assign\n");
    fprintf(out[0],"element id assign\n");

    /* write scalar data files */
    sprintf(name,"ensight/data/engold.%05d.den",jedit);
    out[1]=fopen(name,"w");
    fprintf(out[1],"Per_elem scalar values\n");

    sprintf(name,"ensight/data/engold.%05d.pres",jedit);
    out[2]=fopen(name,"w");
    fprintf(out[2],"Per_elem scalar values\n");
    
    sprintf(name,"ensight/data/engold.%05d.ie",jedit); 
    out[3]=fopen(name,"w");
    fprintf(out[3],"Per_elem scalar values\n");
    
    sprintf(name,"ensight/data/engold.%05d.flg",jedit); 
    out[4]=fopen(name,"w");
    fprintf(out[4],"Per_elem scalar values\n");

    /* write vector data files */
    sprintf(name,"ensight/data/engold.%05d.vel",jedit); 
    out[5]=fopen(name,"w");
    fprintf(out[5],"Per_node vector values\n");

    sprintf(name,"ensight/data/engold.%05d.kfix",jedit); 
    out[6]=fopen(name,"w");
    fprintf(out[6],"Per_node vector values\n");

    
    //printf("jedit= %d l= %d\n",jedit, l);

    for (rank=0;rank<m;rank++) {

      sprintf(name,"ensight_%05d/data/engold.%05d.geo",rank,jedit);
      in=fopen(name,"r");

      for (j=1;j<=8;j++) {

	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}
      }

      fscanf(in,"%d",&kkpll);

      px=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      py=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      pz=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));

      pu=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      pv=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      pw=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));

      kfx1=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      kfx2=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));
      kfx3=(double *) malloc((size_t) ((kkpll+1)*sizeof(double)));

      flg=(int *) malloc((size_t) ((kkpll+1)*sizeof(int)));
      pindx=(int *) malloc((size_t) ((kkpll+1)*sizeof(int)));

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	px[i]=xa;

      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	py[i]=xa;

      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	pz[i]=xa;

      }     

      c=(char)fgetc(in);
      for (j=1;j<=1;j++) {
	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}
      }

      fscanf(in,"%d",&kkzln);
  
      cell=(cell_t *) malloc((size_t) ((kkzln+1)*sizeof(cell_t)));
  
  
      for (i=1;i<=kkzln;i++) {

	fscanf(in,"%d%d%d%d%d%d%d%d",cell[i].pn,cell[i].pn+1,cell[i].pn+2,cell[i].pn+3,
	       cell[i].pn+4,cell[i].pn+5,cell[i].pn+6,cell[i].pn+7);
      }

      fclose(in);


      sprintf(name,"ensight_%05d/data/engold.%05d.vel",rank,jedit);
      in=fopen(name,"r");

      for (j=1;j<=4;j++) {

	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}

      }


      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	pu[i]=xa;

      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	pv[i]=xa;

      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	pw[i]=xa;

      }  

      fclose(in);   

      sprintf(name,"ensight_%05d/data/engold.%05d.kfix",rank,jedit);
      in=fopen(name,"r");

      for (j=1;j<=4;j++) {

	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}
      }


      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	kfx1[i]=xa;

      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	kfx2[i]=xa;
	
      }

      for (i=1;i<=kkpll;i++) {

	fscanf(in,"%le",&xa);
	kfx3[i]=xa;

      }  

      fclose(in);   


    
      sprintf(name,"ensight_%05d/data/engold.%05d.den",rank,jedit);
      in=fopen(name,"r");    


      for (j=1;j<=4;j++) {
	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}

      }

      for (i=1;i<=kkzln;i++) {

	fscanf(in,"%le",&xa);
        cell[i].r=xa;
      }

      fclose(in);
  
      sprintf(name,"ensight_%05d/data/engold.%05d.pres",rank,jedit);
      in=fopen(name,"r");    


      for (j=1;j<=4;j++) {
	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}

      }

      for (i=1;i<=kkzln;i++) {

	fscanf(in,"%le",&xa);
        cell[i].p=xa;
      }

      fclose(in);

      sprintf(name,"ensight_%05d/data/engold.%05d.ie",rank,jedit);
      in=fopen(name,"r");    


      for (j=1;j<=4;j++) {
	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}

      }

      for (i=1;i<=kkzln;i++) {

	fscanf(in,"%le",&xa);
        cell[i].ie=xa;
      }

      fclose(in);


      sprintf(name,"ensight_%05d/data/engold.%05d.flg",rank,jedit);
      in=fopen(name,"r");    


      for (j=1;j<=4;j++) {
	i=0;
	while ((c=(char)fgetc(in))!='\n') {
	  i++;
	}
      }
      for (i=1;i<=kkpll;i++) flg[i]=pindx[i]=0;

      kkzln2=0;
      for (i=1;i<=kkzln;i++) {

	fscanf(in,"%le",&xa);
        cell[i].flg=xa;
        if (xa) {
	  kkzln2++;
          for (j=0;j<8;j++) {
	    flg[cell[i].pn[j]]=1;
	  }
	}
      }

      j=0;
      for (i=1;i<=kkpll;i++) {

	if (flg[i]) {
	  j++;
	  pindx[i]=j;
	}
      }
      
      kkpll2=j;

      fclose(in);  


      fprintf(out[0],"part\n");
      fprintf(out[0],"%10d\n",rank+1);
      fprintf(out[0],"This is the description line for the part\n");
      fprintf(out[0],"coordinates\n");
      fprintf(out[0],"%10lu\n",kkpll2);

      for (i=1;i<=kkpll;i++) 
	if (flg[i]) fprintf(out[0],"%12.5e\n",px[i]);

      for (i=1;i<=kkpll;i++) 
	if (flg[i]) fprintf(out[0],"%12.5e\n",py[i]);

      for (i=1;i<=kkpll;i++) 
	if (flg[i]) fprintf(out[0],"%12.5e\n",pz[i]);

      fprintf(out[0],"hexa8\n");
      fprintf(out[0],"%10lu\n",kkzln2);
  
      for (i=1;i<=kkzln;i++) {
	if (cell[i].flg) {
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[0]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[1]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[2]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[3]]);

	  fprintf(out[0],"%10lu",pindx[cell[i].pn[4]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[5]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[6]]);
	  fprintf(out[0],"%10lu",pindx[cell[i].pn[7]]);
	  fprintf(out[0],"\n");
	}
      }


      for (i=1;i<=4;i++) {
	fprintf(out[i],"part\n");
	fprintf(out[i],"%10d\n",rank+1);
	fprintf(out[i],"hexa8\n");
      }

      for (i=5;i<=6;i++) {
	fprintf(out[i],"part\n");
	fprintf(out[i],"%10d\n",rank+1);
	fprintf(out[i],"coordinates\n");
      }

  
      //write scalar data
      for (i=1;i<=kkzln;i++) {
	if (cell[i].flg) {
	  fprintf(out[1],"%12.5e\n",cell[i].r);
	  fprintf(out[2],"%12.5e\n",cell[i].p);
	  fprintf(out[3],"%12.5e\n",cell[i].ie);
	  fprintf(out[4],"%12.5e\n",cell[i].flg); 
	}// if flag
      }//cell loop

      //write vector data
      for (i=1;i<=kkpll;i++) { 
	if (flg[i]) {
	  fprintf(out[5],"%12.5e\n",pu[i]);
	  fprintf(out[6],"%12.5e\n",kfx1[i]);
	}
      }
      for (i=1;i<=kkpll;i++) { 
	if (flg[i]) {
	  fprintf(out[5],"%12.5e\n",pv[i]);
	  fprintf(out[6],"%12.5e\n",kfx2[i]);
	}
      }
      for (i=1;i<=kkpll;i++) { 
	if (flg[i]) {
	  fprintf(out[5],"%12.5e\n",pw[i]);
	  fprintf(out[6],"%12.5e\n",kfx3[i]);
	}
      }

      printf("done with rank= %d\n",rank+1);

      free(cell);

      free(px);
      free(py);
      free(pz);

      free(pu);
      free(pv);
      free(pw);

      free(kfx1);
      free(kfx2);
      free(kfx3);
      
      free(flg);
      free(pindx);

    } //end rank

    for (i=0;i<=6;i++)
      fclose(out[i]);

    printf("done with jedit= %d\n",jedit);

  } //end jedit


  //write case file
  sprintf(name,"ensight/engold.case");
  out[0]=fopen(name,"w");

  fprintf(out[0],"FORMAT\n");
  fprintf(out[0],"type: ensight gold\n");
  fprintf(out[0],"GEOMETRY\n");
  fprintf(out[0],"model: data/engold.*****.geo\n");
  fprintf(out[0],"VARIABLE\n");
  fprintf(out[0],"scalar per element: den data/engold.*****.den\n");
  fprintf(out[0],"scalar per element: pres data/engold.*****.pres\n");
  fprintf(out[0],"scalar per element: ie data/engold.*****.ie\n");
  fprintf(out[0],"vector per node: vel data/engold.*****.vel\n");
  fprintf(out[0],"vector per node: kfix data/engold.*****.kfix\n");


  fprintf(out[0],"TIME\n");
  fprintf(out[0],"time set: 1\n");
  fprintf(out[0],"number of steps: %4d\n",n+1);
  fprintf(out[0],"filename numbers: \n");
  for (i=0;i<=n;i++)
    fprintf(out[0],"%05d\n",i);

  fprintf(out[0],"time values: \n");
  for (i=0;i<=n;i++)
    fprintf(out[0],"%12.5e\n",i*etime);

  fclose(out[0]);


}

