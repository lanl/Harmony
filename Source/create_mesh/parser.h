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
typedef struct _mat {
  double cv;
  double g;
  double csmin;
  double pmin;
} mat_t;

typedef struct _contour {
  char name[50];
  int points;
} contour_t;

typedef struct _contour_data {
  contour_t contour[20];
  int NC;
} con_data_t;

typedef struct _fill {
  int type;
  int j1,j2;
  int i1,i2;
  int k1,k2;
  double x1,x2;
  double y1,y2;
  double z1,z2;
  int con1,con2;
  double rad1,rad2;
  double theta1,theta2;
  double phi1,phi2;
  int reg1,reg2;
  double u,v,w;
  double en,rho;
} fill_t;

typedef struct _fill_data {
  fill_t fill[40];
  int NF;
} fill_data_t;


typedef struct _bound {
  int ifix[3];
  int planex,planey,planez;
  int cylinder;
  int uniform;
  double x,y,z;
} bound_t;

typedef struct _bound_data {

  int NB;
  bound_t bnd[40];

} bound_data_t;

typedef struct _rlx_data { //data for mesh relaxation specifications

  int type;
  int nlag;         //number of logical j lines that are specified to be Lagrangian
  int jlines[15];   //list of logical j lines that are specified to be Lagrangian

} rlx_data_t;


typedef struct _global_data {
  int iorder;
  int ihour;
  int edit_type;
  int iremap;
  int isource;
  double TFINAL;
  double DTMIN;
  double DT;
  int igeo;
  double freeze_time;
  double relax_time;
  double etime;
  double CA;
  double C1;
  double C2;
  int LIMX;
  int LIMY;
  int LIMZ;
  int NPEX;
  int NPEY;
  int NPEZ;    
  int max_steps;
  int input_type;
  int kefix;
  int ihgforce;
    
  double time;
  int ifstop;
  int jedit;
  int ncyc;
  int trilink_spin;    
} global_data_t;

// make material data structure
mat_t mat;

// make contour data structure
con_data_t cd;

// make fill instruction data structure
fill_data_t fd;

// make global data structure
global_data_t gd;

// make boundary instruction data structure
bound_data_t bd;

//ALE mesh relaxation data structure
rlx_data_t rlx;

int parser(char *file)
{

  FILE *in;

  char c,buffer[100000],word[1000];

  int ii,i,j,l,k,jj;
  char *string[10000];
  double x,y;

  int ng=0;


  /* read the input file */
  in=fopen(file, "r");

  i=0;
  while ((buffer[i]=(char)fgetc(in))!=EOF)
    i++;

  buffer[i]='\0';

  fclose(in);


  /* determine the number of strings, k,  and store each string */
  i=0;
  k=0;

  c=buffer[i];

  while (c!='\0') {

    while (c==' ' || c=='\t' || c=='=' || c=='\n' || c=='\r') {
      i++;
      c=buffer[i];
    }
    if (c!='\0') {
      l=0;
      word[l]=c;
      i++;
      c=buffer[i];

      while (c!=' ' && c!='\t' && c!='=' && c!='\n' && c!='\r') {
	l++;
	word[l]=c;
	i++;
	c=buffer[i];
      }
      l++;
      word[l]='\0';

      string[k]=(char *)malloc((size_t)((l+1)*sizeof(char)));
      for (j=0;j<=l;j++)
	string[k][j]=word[j];
      k++;
    }

  }
  k--;

  /* find the number of fill instructions */

  l=0;
  for (j=0;j<=k;j++)
    if (strcmp("fill",string[j])==0) l++;

  printf("Number of fill instructions=%d\n",l);
  fd.NF=l;

  /* find the number of contours */

  l=0;
  for (j=0;j<=k;j++)
    if (strcmp("contour",string[j])==0) l++;

  printf("Number of contours=%d\n",l);

  cd.NC=l;

  /* find the number of boundary specification */

  l=0;
  for (j=0;j<=k;j++)
    if (strcmp("bound",string[j])==0) l++;

  printf("Number of boundary specifications=%d\n",l);

  bd.NB=l;



  j=0;

  while (j<k) {

    l=j;

    if (strcmp("global",string[j])==0) {
      j++;
      while (strcmp("end",string[j])!=0) {
	if (strcmp("iorder",string[j])==0) {
	  j++;
	  gd.iorder=atoi(string[j]);
	  printf("iorder=%d\n",gd.iorder);
	}
	else if (strcmp("ihour",string[j])==0) {
	  j++;
	  gd.ihour=atoi(string[j]);
	  printf("hourglass flag=%d\n",gd.ihour);
	}
	else if (strcmp("edit_type",string[j])==0) {
	  j++;
	  gd.edit_type=atoi(string[j]);
	  printf("edit type=%d\n",gd.edit_type);
	}
	else if (strcmp("iremap",string[j])==0) {
	  j++;
	  gd.iremap=atoi(string[j]);
	  printf("remap flag=%d\n",gd.iremap);
	}
	else if (strcmp("isource",string[j])==0) {
	  j++;
	  gd.isource=atoi(string[j]);
	  printf("isource flag=%d\n",gd.isource);
	}
	else if (strcmp("tstop",string[j])==0) {
	  j++;
	  gd.TFINAL=atof(string[j]);
	  printf("tstop=%e\n",gd.TFINAL);
	}
	else if (strcmp("dtmin",string[j])==0) {
	  j++;
	  gd.DTMIN=atof(string[j]);
	  printf("dtmin=%e\n",gd.DTMIN);
	}
	else if (strcmp("dt",string[j])==0) {
	  j++;
	  gd.DT=atof(string[j]);
	  printf("dt=%e\n",gd.DT);
	}
	else if (strcmp("igeo",string[j])==0) {
	  j++;
	  gd.igeo=atoi(string[j]);
	  printf("geometry option=%d\n",gd.igeo);
	}
	else if (strcmp("freeze",string[j])==0) {
	  j++;
	  gd.freeze_time=atof(string[j]);
	  printf("Mesh freeze time=%e\n",gd.freeze_time);
	}
	else if (strcmp("relax",string[j])==0) {
	  j++;
	  gd.relax_time=atof(string[j]);
	  printf("Mesh relax time=%e\n",gd.relax_time);
	}
	else if (strcmp("etime",string[j])==0) {
	  j++;
	  gd.etime=atof(string[j]);
	  printf("Edit time interval=%e\n",gd.etime);
	}
	else if (strcmp("CA",string[j])==0) {
	  j++;
	  gd.CA=atof(string[j]);
	  printf("Hourglass viscosity coefficient=%e\n",gd.CA);
	}
	else if (strcmp("C1",string[j])==0) {
	  j++;
	  gd.C1=atof(string[j]);
	  printf("Quadratic viscosity coefficient=%e\n",gd.C1);
	}
	else if (strcmp("C2",string[j])==0) {
	  j++;
	  gd.C2=atof(string[j]);
	  printf("Liner viscosity coefficient=%e\n",gd.C2);
	}
	else if (strcmp("LIMX",string[j])==0) {
	  j++;
	  gd.LIMX=atoi(string[j]);
	  printf("Number of x zones=%d\n",gd.LIMX);
	}
	else if (strcmp("LIMY",string[j])==0) {
	  j++;
	  gd.LIMY=atoi(string[j]);
	  printf("Number of y zones=%d\n",gd.LIMY);
	}
	else if (strcmp("LIMZ",string[j])==0) {
	  j++;
	  gd.LIMZ=atoi(string[j]);
	  printf("Number of z zones=%d\n",gd.LIMZ);
	}
	else if (strcmp("NPEX",string[j])==0) {
	  j++;
	  gd.NPEX=atoi(string[j]);
	  printf("Number of processors in x=%d\n",gd.NPEX);
	}
	else if (strcmp("NPEY",string[j])==0) {
	  j++;
	  gd.NPEY=atoi(string[j]);
	  printf("Number of processors in y=%d\n",gd.NPEY);
	}
	else if (strcmp("NPEZ",string[j])==0) {
	  j++;
	  gd.NPEZ=atoi(string[j]);
	  printf("Number of processors in z=%d\n",gd.NPEZ);
	}	
	else if (strcmp("max_steps",string[j])==0) {
	  j++;
	  gd.max_steps=atoi(string[j]);
	  printf("Maximum number of cycles = %d\n",gd.max_steps);
	}
	else if (strcmp("input_type",string[j])==0) {
	  j++;
	  gd.input_type=atoi(string[j]);
	  printf("Input Type = %d\n",gd.input_type);
	}
	else if (strcmp("kefix",string[j])==0) {
	  j++;
	  gd.kefix=atoi(string[j]);
	  printf("kefix flag = %d\n",gd.kefix);
	}
	else if (strcmp("ihgforce",string[j])==0) {
	  j++;
	  gd.ihgforce=atoi(string[j]);
	  printf("hourglass flag = %d\n",gd.ihgforce);
	}
	j++;
      }
      if (j!=k) j++;
    }

    if (strcmp("relax",string[j])==0) {
      printf("ALE mesh relaxation specified\n");
      j++;

      while (strcmp("end",string[j])!=0) {
	if (strcmp("type",string[j])==0) {
	  j++;
	  rlx.type=atoi(string[j]);
	  printf("ALE relaxer type =%d\n",rlx.type);
	}
	else if (strcmp("nlag",string[j])==0) {
	  j++;
	  rlx.nlag=atoi(string[j]);
	  printf("Number of Lagragian interfaces =%d\n",rlx.nlag);
	}
	else if (strcmp("J_lines",string[j])==0) {
	  
	  for (ii=0;ii<rlx.nlag;ii++) {
	    j++;
	    rlx.jlines[ii]=atoi(string[j]);
	  }
	  
	  printf("Lagrangian Interfaces on J lines ");
	  for (ii=0;ii<rlx.nlag;ii++) {
	    printf("%d ",rlx.jlines[ii]);
	  }
	  printf("\n");
	}

	j++;
      }
      if (j!=k) j++;
    }



    if (strcmp("bound",string[j])==0) {
      j++;
      i=atoi(string[j])-1;
      bd.bnd[i].planex=bd.bnd[i].planey=bd.bnd[i].planez=0;
      bd.bnd[i].cylinder=0;
      bd.bnd[i].ifix[0]=bd.bnd[i].ifix[1]=bd.bnd[i].ifix[2];
      printf("Boundary specification number=%d\n",i);
      j++;

      while (strcmp("end",string[j])!=0) {
	if (strcmp("planex",string[j])==0) {
	  j++;
	  bd.bnd[i].planex=atoi(string[j]);
	  printf("x plane condition =%d\n",bd.bnd[i].planex);
	}
	else if (strcmp("planey",string[j])==0) {
	  j++;
	  bd.bnd[i].planey=atoi(string[j]);
	  printf("y plane condition =%d\n",bd.bnd[i].planey);
	}
	else if (strcmp("planez",string[j])==0) {
	  j++;
	  bd.bnd[i].planez=atoi(string[j]);
	  printf("z plane condition =%d\n",bd.bnd[i].planez);
	}
	else if (strcmp("cylinder",string[j])==0) {
	  j++;
	  bd.bnd[i].cylinder=atoi(string[j]);
	  printf("cylinder condition =%d\n",bd.bnd[i].cylinder);
	}
	else if (strcmp("uniform",string[j])==0) {
	  j++;
	  bd.bnd[i].uniform=atoi(string[j]);
	  printf("uniform condition =%d\n",bd.bnd[i].uniform);
	}	
	else if (strcmp("x",string[j])==0) {
	  j++;
	  bd.bnd[i].x=atof(string[j]);
	  printf("x=%e\n",bd.bnd[i].x);
	}
	else if (strcmp("y",string[j])==0) {
	  j++;
	  bd.bnd[i].y=atof(string[j]);
	  printf("y=%e\n",bd.bnd[i].y);
	}
	else if (strcmp("z",string[j])==0) {
	  j++;
	  bd.bnd[i].z=atof(string[j]);
	  printf("z=%e\n",bd.bnd[i].z);
	}
	else if (strcmp("ifix1",string[j])==0) {
	  j++;
	  bd.bnd[i].ifix[0]=atoi(string[j]);
	  printf("ifix(1)=%d\n",bd.bnd[i].ifix[0]);
	}
	else if (strcmp("ifix2",string[j])==0) {
	  j++;
	  bd.bnd[i].ifix[1]=atoi(string[j]);
	  printf("ifix(2)=%d\n",bd.bnd[i].ifix[1]);
	}
	else if (strcmp("ifix3",string[j])==0) {
	  j++;
	  bd.bnd[i].ifix[2]=atoi(string[j]);
	  printf("ifix(3)=%d\n",bd.bnd[i].ifix[2]);
	}

	j++;
      }
      if (j!=k) j++;
    }


    if (strcmp("fill",string[j])==0) {
      j++;
      i=atoi(string[j])-1;
      printf("Fill instruction number=%d\n",i);
      fd.fill[i].rad1=fd.fill[i].rad2=0.0;
      j++;

      while (strcmp("end",string[j])!=0) {
	if (strcmp("type",string[j])==0) {
	  j++;
	  fd.fill[i].type=atoi(string[j]);
	  printf("Type=%d\n",fd.fill[i].type);
	}
	else if (strcmp("contour1",string[j])==0) {
	  j++;
	  fd.fill[i].con1=atoi(string[j])-1;
	  printf("Contour 1=%d\n",fd.fill[i].con1);
	}
	else if (strcmp("contour2",string[j])==0) {
	  j++;
	  fd.fill[i].con2=atoi(string[j])-1;
	  printf("Contour 2=%d\n",fd.fill[i].con2);
	}
	else if (strcmp("j1",string[j])==0) {
	  j++;
	  fd.fill[i].j1=atoi(string[j]);
	  printf("j1=%d\n",fd.fill[i].j1);
	}
	else if (strcmp("j2",string[j])==0) {
	  j++;
	  fd.fill[i].j2=atoi(string[j]);
	  printf("j2=%d\n",fd.fill[i].j2);
	}
	else if (strcmp("k1",string[j])==0) {
	  j++;
	  fd.fill[i].k1=atoi(string[j]);
	  printf("k1=%d\n",fd.fill[i].k1);
	}
	else if (strcmp("k2",string[j])==0) {
	  j++;
	  fd.fill[i].k2=atoi(string[j]);
	  printf("k2=%d\n",fd.fill[i].k2);
	}
	else if	(strcmp("i1",string[j])==0) {
	  j++;
	  fd.fill[i].i1=atoi(string[j]);
	  printf("i1=%d\n",fd.fill[i].i1);
	}
	else if (strcmp("i2",string[j])==0) {
	  j++;
	  fd.fill[i].i2=atoi(string[j]);
	  printf("i2=%d\n",fd.fill[i].i2);
	}
	else if (strcmp("region1",string[j])==0) {
	  j++;
	  fd.fill[i].reg1=atoi(string[j])-1;
	  printf("Background region=%d\n",fd.fill[i].reg1);
	}
	else if (strcmp("region2",string[j])==0) {
	  j++;
	  fd.fill[i].reg2=atoi(string[j])-1;
	  printf("Painted region=%d\n",fd.fill[i].reg2);
	}
	else if (strcmp("theta1",string[j])==0) {
	  j++;
	  fd.fill[i].theta1=atof(string[j]);
	  printf("Theta 1=%e\n",fd.fill[i].theta1);
	}
	else if (strcmp("theta2",string[j])==0) {
	  j++;
	  fd.fill[i].theta2=atof(string[j]);
	  printf("Theta 2=%e\n",fd.fill[i].theta2);
	}
	else if (strcmp("phi1",string[j])==0) {
	  j++;
	  fd.fill[i].phi1=atof(string[j]);
	  printf("Phi 1=%e\n",fd.fill[i].phi1);
	}
	else if (strcmp("phi2",string[j])==0) {
	  j++;
	  fd.fill[i].phi2=atof(string[j]);
	  printf("Phi 2=%e\n",fd.fill[i].phi2);
	}
	else if (strcmp("radius1",string[j])==0) {
	  j++;
	  fd.fill[i].rad1=atof(string[j]);
	  printf("Radius 1=%e\n",fd.fill[i].rad1);
	}
	else if (strcmp("radius2",string[j])==0) {
	  j++;
	  fd.fill[i].rad2=atof(string[j]);
	  printf("Radius 2=%e\n",fd.fill[i].rad2);
	}
	else if	(strcmp("x1",string[j])==0) {
	  j++;
	  fd.fill[i].x1=atof(string[j]);
	  printf("x1=%e\n",fd.fill[i].x1);
	}
	else if	(strcmp("x2",string[j])==0) {
	  j++;
	  fd.fill[i].x2=atof(string[j]);
	  printf("x2=%e\n",fd.fill[i].x2);
	}
	else if	(strcmp("y1",string[j])==0) {
	  j++;
	  fd.fill[i].y1=atof(string[j]);
	  printf("y1=%e\n",fd.fill[i].y1);
	}
	else if (strcmp("y2",string[j])==0) {
	  j++;
	  fd.fill[i].y2=atof(string[j]);
	  printf("y2=%e\n",fd.fill[i].y2);
	}
	else if	(strcmp("z1",string[j])==0) {
	  j++;
	  fd.fill[i].z1=atof(string[j]);
	  printf("z1=%e\n",fd.fill[i].z1);
	}
	else if (strcmp("z2",string[j])==0) {
	  j++;
	  fd.fill[i].z2=atof(string[j]);
	  printf("z2=%e\n",fd.fill[i].z2);
	}
	else if (strcmp("energy",string[j])==0) {
	  j++;
	  fd.fill[i].en=atof(string[j]);
	  printf("energy=%e\n",fd.fill[i].en);
	}
	else if	(strcmp("rho",string[j])==0) {
	  j++;
	  fd.fill[i].rho=atof(string[j]);
	  printf("rho=%e\n",fd.fill[i].rho);
	}
	else if (strcmp("u",string[j])==0) {
	  j++;
	  fd.fill[i].u=atof(string[j]);
	  printf("u=%e\n",fd.fill[i].u);
	}
	else if (strcmp("v",string[j])==0) {
	  j++;
	  fd.fill[i].v=atof(string[j]);
	  printf("v=%e\n",fd.fill[i].v);
	}
	else if (strcmp("w",string[j])==0) {
	  j++;
	  fd.fill[i].w=atof(string[j]);
	  printf("w=%e\n",fd.fill[i].w);
	}	
	j++;
      }
      if (j!=k) j++;
    }


    if (strcmp("contour",string[j])==0) {
      j++;
      i=atoi(string[j])-1;
      printf("Contour number=%d\n",i);
      j++;
      while (strcmp("end",string[j])!=0) {

	if (strcmp("file",string[j])==0) {
	  j++;
	  ii=0;
	  while (string[j][ii]!='\0') {
	    cd.contour[i].name[ii]=string[j][ii];
	    ii++;
	  }
	  cd.contour[i].name[ii]='\0';
	  printf("Contour File Name=%s\n",cd.contour[i].name);
	}
	else if (strcmp("points",string[j])==0) {
	  j++;
	  cd.contour[i].points=atof(string[j]);
	  printf("Number of points on the contour=%d\n",cd.contour[i].points);
	  cd.contour[i].points--;
	}
	j++;
      }
      if (j!=k) j++;
    }


    if (strcmp("eos",string[j])==0) {
      printf("%s\n",string[j]);
      j++;
      while (strcmp("end",string[j])!=0) {
	if (strcmp("cv",string[j])==0) {
	  j++;
	  mat.cv=atof(string[j]);
	  printf("cv=%e\n",mat.cv);
	}
	else if	(strcmp("gam",string[j])==0) {
	  j++;
	  mat.g=atof(string[j]);
	  printf("gam=%e\n",mat.g);
	}
	else if (strcmp("csmin",string[j])==0) {
	  j++;
	  mat.csmin=atof(string[j]);
	  printf("Minimum Sound Speed=%e\n",mat.csmin);
	  mat.csmin*=mat.csmin;
	}
	else if (strcmp("pmin",string[j])==0) {
	  j++;
	  mat.pmin=atof(string[j]);
	  printf("Minimum Pressure=%e\n",mat.pmin);
	}	
	j++;
      }


      if (j!=k) j++;
    }

    if (l==j) {
      printf("%s is not understood\n" ,string[j]);
      gd.ifstop=1;
      break;
    }
  }


  for (j=0;j<=k;j++)
    free(string[j]);

  return 0;

}
