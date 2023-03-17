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
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <array>

class my_cell_t {
  
public :
  int pn[8];
  double xc,yc,zc;
  double flg,rank;
  int proc,pos;
  int num;
  int cid;
  int mype;

  my_cell_t();
};


my_cell_t::my_cell_t()
{
  int i;
  
  xc=yc=zc=0.0;
  flg=rank=0.0;
  num=0;
  cid=-1;
  mype=-1;
  for (i=0;i<8;i++)
    pn[i]=-1;
  proc=pos=-1;

}


class my_vertex_t {
  
public :
  double px[3];

  int flg;

  my_vertex_t();
};

my_vertex_t::my_vertex_t()
{
  int i;
  
  for (i=0;i<3;i++) {
    px[i]=0.0;
  }

  flg=0;

}



int main(int argc, char** argv)
{

  int i,j,k,l,id,m,n,p,kk,flg2;

  int ifstr,jmat,lmat;

  int LIMX,LIMY,LIMZ,JT;

  int kkzln,kkzln2,kkpll,kkpll2;

  FILE *in,*in2,*out[100];

  std::vector<my_cell_t> cell;
  std::vector<my_vertex_t> vertex;
  
  std::vector<my_cell_t> cell0;
  std::vector<my_vertex_t> vertex0;

  char name[100],c;

  double pres,xa,ya,za,ua,va,wa;

  double x1,x2,x3,x4,y1,y2,y3,y4;
  double z1,z2,z3,z4;

  double sum,sum1;

  double Pi=2.*acos(0.);

  int jedit,rank;

  double edtimes[250];

  int ifcch, *flag2;

  int num1, *send_rank, *send_cell, *recv_rank, *recv_cell;

  int i1,i2,j1,j2,k1,k2;
  double dx,dy,dz;
  double theta1,rad;
  double xc1[3];
 
  int count,mpi_count;

  int idebug=0;

  m=atoi(argv[1]);

  jedit=0;
  //open master mesh file
  sprintf(name,"mesh.geo");
  in=fopen(name,"r");
  for (j=1;j<=12;j++) {

    i=0;
    while ((c=(char)fgetc(in))!='\n') {
        i++;
    }
  }

  fscanf(in,"%d",&kkpll);
  vertex0.reserve(kkpll);
  vertex0.resize(kkpll);


  //printf("size of vertex array= %d %d\n",kkpll,vertex0.size());
  
  for (i=0;i<kkpll;i++) {

    fscanf(in,"%d",&k);

  }
    
  for (i=0;i<kkpll;i++) {

    fscanf(in,"%le",&xa);
    vertex0[i].px[0]=xa;

  }

  for (i=0;i<kkpll;i++) {

    fscanf(in,"%le",&xa);
    vertex0[i].px[1]=xa;

  }

  for (i=0;i<kkpll;i++) {

    fscanf(in,"%le",&xa);
    vertex0[i].px[2]=xa;

  }


  c=(char)fgetc(in);
  for (j=1;j<=1;j++) {
    i=0;
    while ((c=(char)fgetc(in))!='\n') {
      i++;
    }
  }

  fscanf(in,"%d",&kkzln);
  cell0.reserve(kkzln);
  cell0.resize(kkzln);
  

  //printf("size of cell array= %d %d\n",kkzln,cell0.size());

  for (i=0;i<kkzln;i++) {
      
    fscanf(in,"%d",&k);

  } 
  
  for (i=0;i<kkzln;i++) {

    fscanf(in,"%d%d%d%d%d%d%d%d",cell0[i].pn,cell0[i].pn+1,cell0[i].pn+2,cell0[i].pn+3,
	   cell0[i].pn+4,cell0[i].pn+5,cell0[i].pn+6,cell0[i].pn+7);

    for (j=0;j<8;j++)
      cell0[i].pn[j]--;

  }

  fclose(in);

  //done reading global mesh
  
  count=0; //count total number of ghost cells over all ranks
  for (rank=0;rank<m;rank++) { //loop over all ranks, read mesh for each rank

    sprintf(name,"ensight_%05d/engold.geo",rank);

    in=fopen(name,"r");

    for (j=1;j<=12;j++) {

      i=0;
      while ((c=(char)fgetc(in))!='\n') {
	i++;
      }

    }

    fscanf(in,"%d",&kkpll);
    vertex.clear();
    vertex.reserve(kkpll);
    vertex.resize(kkpll);


    //printf("size of vertex array= %d %d\n",kkpll,vertex.size());

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%d",&k);

    }   

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[0]=xa;

    }

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[1]=xa;

    }

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[2]=xa;

    }


    c=(char)fgetc(in);
    for (j=1;j<=1;j++) {
      i=0;
      while ((c=(char)fgetc(in))!='\n') {
        i++;
      }
    }

    fscanf(in,"%d",&kkzln);
    cell.clear();
    cell.reserve(kkzln);
    cell.resize(kkzln);


    //printf("size of cell array= %d %d\n",kkzln,cell.size());

    for (i=0;i<kkzln;i++) {
      
      fscanf(in,"%d",&k);

    }    
  
    for (i=0;i<kkzln;i++) {

      fscanf(in,"%d%d%d%d%d%d%d%d",cell[i].pn,cell[i].pn+1,cell[i].pn+2,cell[i].pn+3,
	     cell[i].pn+4,cell[i].pn+5,cell[i].pn+6,cell[i].pn+7);

      for (j=0;j<8;j++)
	cell[i].pn[j]--;      
      
    }

    fclose(in);

    //done reading geo file



    //read ownership flag
    sprintf(name,"ensight_%05d/engold.flg",rank);
    //printf("%s\n",name);
    in=fopen(name,"r");   

    for (i=0;i<kkzln;i++) {

      fscanf(in,"%d",&k);
      cell[i].flg=1.0*k;

      if (cell[i].flg<0.9) count++; //count ghost cell
 
    }

    fclose(in);

    //read global id for each cell
    sprintf(name,"ensight_%05d/engold.cid",rank);
    //printf("%s\n",name);
    in=fopen(name,"r");

    for (int c=0;c<kkzln;c++) {

      fscanf(in,"%d",&k);
      cell[c].cid=k;

      //store rank and local cell index c on the master mesh
      if (cell[c].flg>0.9) {
	cell0[k].mype=rank;
	cell0[k].pos=c;
      }    
      
    }
    fclose(in);
    
  
    //printf("done with rank= %d\n",rank+1);

  } //end first loop over ranks

  //write the rank that owns each cell on the global mesh
  sprintf(name,"mesh.rank");
  in=fopen(name,"w");

  fprintf(in,"Per_elem scalar values\n");
  fprintf(in,"part\n");
  fprintf(in,"%10d\n",1);
  fprintf(in,"hexa8\n");

  for (i=0;i<cell0.size();i++) {
      
    fprintf(in,"%12.5e\n",1.0*cell0[i].mype);
  }
  fclose(in);
  
  //update the case file for the global mesh
  sprintf(name,"mesh.case");
  in=fopen(name,"w");
  fprintf(in,"# Wed Aug 29 23:50:44 2018\n");
  fprintf(in,"# EnSight Gold Model\n");
  fprintf(in,"# Produced with EnSight Ouput API, version 1.0.4\n");
  fprintf(in,"\n");
  fprintf(in,"FORMAT\n");
  fprintf(in,"type:                   ensight gold\n");
  fprintf(in,"\n");
  fprintf(in,"GEOMETRY\n");
  fprintf(in,"model:                  mesh.geo\n");
  fprintf(in,"VARIABLE\n");
  fprintf(in,"scalar per element: rank mesh.rank\n");
  fclose(in);  

  //create a list of all MPI cell to cell communications, including rank and local cell id for both send and receive cells

  recv_cell= new int[count];
  send_cell= new int[count];
  recv_rank= new int[count];
  send_rank= new int[count];

  mpi_count=count;

    
  //second loop over ranks
  count=0;
  for (rank=0;rank<m;rank++) {

    sprintf(name,"ensight_%05d/engold.geo",rank);
    in=fopen(name,"r");

    for (j=1;j<=12;j++) {

      i=0;
      while ((c=(char)fgetc(in))!='\n') {
	i++;
      }

    }

    fscanf(in,"%d",&kkpll);
    vertex.clear();
    vertex.reserve(kkpll);
    vertex.resize(kkpll);


    //printf("size of vertex array= %d %d\n",kkpll,vertex.size());

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%d",&k);

    }    

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[0]=xa;

    }

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[1]=xa;

    }

    for (i=0;i<kkpll;i++) {

      fscanf(in,"%le",&xa);
      vertex[i].px[2]=xa;

    }


    c=(char)fgetc(in);
    for (j=1;j<=1;j++) {
      i=0;
      while ((c=(char)fgetc(in))!='\n') {
        i++;
      }
    }

    fscanf(in,"%d",&kkzln);
    cell.clear();
    cell.reserve(kkzln);
    cell.resize(kkzln);

    //printf("size of cell array= %d %d\n",kkzln,cell.size());

    for (i=0;i<kkzln;i++) {
      
      fscanf(in,"%d",&k);

    }      
  
    for (i=0;i<kkzln;i++) {

      fscanf(in,"%d%d%d%d%d%d%d%d",cell[i].pn,cell[i].pn+1,cell[i].pn+2,cell[i].pn+3,
	     cell[i].pn+4,cell[i].pn+5,cell[i].pn+6,cell[i].pn+7);

      for (j=0;j<8;j++)
	cell[i].pn[j]--;      
      
    }

    fclose(in);

    //done reading geo file

    //read ownership flag
    sprintf(name,"ensight_%05d/engold.flg",rank);
    //printf("%s\n",name);
    in=fopen(name,"r");   

    for (i=0;i<kkzln;i++) {

      fscanf(in,"%d",&k);
      cell[i].flg=1.0*k;

    }

    fclose(in);

    //read global id for each cell
    sprintf(name,"ensight_%05d/engold.cid",rank);
    //printf("%s\n",name);
    in=fopen(name,"r");

    for (int c=0;c<kkzln;c++) {

      fscanf(in,"%d",&k);
      cell[c].cid=k;
      
    }
    fclose(in);

    kkzln2=cell0.size();

    //loop over all cells to identify ghost cells and associated MPI cell to cell communications
    for (i=0;i<kkzln;i++) {

      cell[i].rank=1.0*rank;

      if (cell[i].flg<0.9) { //ghost cell
    

	recv_cell[count]=i;
	recv_rank[count]=rank;
	
	k=cell[i].cid;
	send_rank[count]=cell0[k].mype;
	send_cell[count]=cell0[k].pos;
	cell[i].rank=1.0*cell0[k].mype;
	count++;
	    
      } //end if ghost cell

    } //end cell loop


    //write cell rank
    sprintf(name,"ensight_%05d/engold.rank",rank);
    in=fopen(name,"w"); 

    fprintf(in,"Per_elem scalar values\n");
    fprintf(in,"part\n");
    fprintf(in,"%10d\n",1);
    fprintf(in,"hexa8\n");

    for (i=0;i<kkzln;i++) {
      
      fprintf(in,"%12.5e\n",cell[i].rank);
    }

    fclose(in);

    //write cell ownership flag
    sprintf(name,"ensight_%05d/engold.flag",rank);
    in=fopen(name,"w"); 

    fprintf(in,"Per_elem scalar values\n");
    fprintf(in,"part\n");
    fprintf(in,"%10d\n",1);
    fprintf(in,"hexa8\n");

    for (i=0;i<kkzln;i++) {
      
      fprintf(in,"%12.5e\n",cell[i].flg);
    }

    fclose(in);

    //update the case file for each rank
    sprintf(name,"ensight_%05d/engold.case",rank);
    in=fopen(name,"w");
    fprintf(in,"# Wed Aug 29 23:50:44 2018\n");
    fprintf(in,"# EnSight Gold Model\n");
    fprintf(in,"# Produced with EnSight Ouput API, version 1.0.4\n");
    fprintf(in,"\n");
    fprintf(in,"FORMAT\n");
    fprintf(in,"type:                   ensight gold\n");
    fprintf(in,"\n");
    fprintf(in,"GEOMETRY\n");
    fprintf(in,"model:                  engold.geo\n");
    fprintf(in,"VARIABLE\n");
    fprintf(in,"scalar per element: rank engold.rank\n");
    fprintf(in,"scalar per element: flag engold.flag\n");
    fclose(in);   
    
  } //end second loop over ranks

  //third loop over ranks, write connectivity file for each rank
  for (rank=0;rank<m;rank++) {  

    //print the process communication table
    sprintf(name,"ensight_%05d/connectivity.dat",rank);
    in=fopen(name,"w");

    if (idebug) sprintf(name,"ensight_%05d/connectivity_debug.dat",rank);
    if (idebug) in2=fopen(name,"w");   

    //print list of cells sent by rank, ordered by rank
    for (i=0;i<m;i++) {
      
      //find number of cells send by rank to process i
      int local_count=0;
      for (k=0;k<count;k++) //loop over all MPI cell to cell communications
	if (send_rank[k]==rank && recv_rank[k]==i) local_count++;

      fprintf(in,"%d\n",local_count);
      if (idebug) fprintf(in2,"number of cells sent to process %d equals %d\n",i,local_count);
      
      for (k=0;k<count;k++) //loop over all MPI cell to cell communications
	if (send_rank[k]==rank && recv_rank[k]==i) {
	  fprintf(in,"%d %d\n",send_cell[k],recv_cell[k]);
	  if (idebug) fprintf(in2,"%d %d\n",send_cell[k],recv_cell[k]);
	}
      
    }        

    //print list of cells received by rank, ordered by rank
    for (i=0;i<m;i++) {

      //find number of cells received from process i by rank
      int local_count=0;
      for (k=0;k<count;k++) //loop over all MPI cell to cell communications
	if (send_rank[k]==i && recv_rank[k]==rank) local_count++;
      
      fprintf(in,"%d\n",local_count);
      if (idebug) fprintf(in2,"number of cells received from process %d equals %d\n",i,local_count);
      
      for (k=0;k<count;k++) //loop over all MPI cell to cell communications
	if (send_rank[k]==i && recv_rank[k]==rank) {
	  fprintf(in,"%d %d\n",recv_cell[k],send_cell[k]);
	  if (idebug) fprintf(in2,"%d %d\n",recv_cell[k],send_cell[k]);
	}

    }
  
    fclose(in);
    if (idebug) fclose(in2);

  } //end loop over ranks

  delete[] send_rank;
  delete[] recv_rank;
  delete[] send_cell;
  delete[] recv_cell;

}

