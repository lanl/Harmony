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

#include "mesh.h"

int main(int argc, char** argv)
{
  char buff[100]="engold.geo";

  mesh_t mesh(buff);

  FILE *out;

  int i,j,k,n,z,p,m,ii;

  int p1,p2,idup,j1,j2,mlist[4],list8[8];

  int *vflag;

  out=fopen("mesh_file.dat","w");

  fprintf(out,"%10d\n",mesh.kkpll);

  //for each point determine a list of edges
  vflag= new int[6*mesh.kkpll];

  for (i=1;i<=mesh.kkpll;i++)
    for (n=1;n<=6;n++)
      vflag[findex(6,i,n)]=-1;

  for (i=1;i<=mesh.kkell;i++) {

    p1=mesh.kkep1[findex1(i)];

    j=1;
    while (vflag[findex(6,p1,j)]!=-1)
      j=j+1;

    vflag[findex(6,p1,j)]=i;

    p2=mesh.kkep2[findex1(i)];          
          
    j=1;
    while (vflag[findex(6,p2,j)]!=-1)
      j=j+1;
 
    vflag[findex(6,p2,j)]=i;  


  }

  for (j=1;j<=mesh.kkpll;j++){
    fprintf(out,"%10d",j);
    for (n=1;n<=6;n++)
      fprintf(out,"%10d",vflag[findex(6,j,n)]);
    fprintf(out,"\n");
  }

  delete[] vflag;

  //for each point determine a list of faces
  vflag=new int[12*mesh.kkpll];
  
  for (i=1;i<=mesh.kkpll;i++)
    for (n=1;n<=12;n++)
      vflag[findex(12,i,n)]=-1;
  
  for (i=1;i<=mesh.kkfll;i++) {

    //list of points for this face
    for (k=1;k<=4;k++)
      mlist[findex1(k)]=-1;
          
    //cell corresponding to this face
    z=mesh.kkfz[findex1(i)];

    if (z>mesh.kkzln) printf("cell id= %d\n",z);
 
          
    //first side in the cell
    j1=mesh.kkzs[findex1(z)]-23;

    //last side in the cell
    j2=mesh.kkzs[findex1(z)];

    //loop over sides in the cell 
    for (j=j1;j<=j2;j++) {

      if (mesh.kksf[findex1(j)]==i) {//this side belongs to the face

	//assign p1 to the face
	p1=mesh.kksp1[findex1(j)];

	idup=0;

	k=1;
	while (mlist[findex1(k)]!=-1) {

	  if (mlist[findex1(k)]==p1) idup=1;
                      
	  k=k+1;

	  if (k==5) break;
	}
                   
	if (idup==0 && k<=4) mlist[findex1(k)]=p1;
                
	//assign p2 to the face
	p2=mesh.kksp2[findex1(j)];

	idup=0;

	k=1;
	while (mlist[findex1(k)]!=-1) {

	  if (mlist[findex1(k)]==p2) idup=1;
                      
	  k=k+1;

	  if (k==5) break;
	}	

	if (idup==0 && k<=4) mlist[findex1(k)]=p2;

      } //end if                     
    }

    //tag all points on the face
    for (k=1;k<=4;k++) {

      p1=mlist[findex1(k)];
      j=1;
      while (vflag[findex(12,p1,j)]!=-1)
	j=j+1;   
              
      vflag[findex(12,p1,j)]=i;

      if (p1==-1) printf("p1= -1 on face\n");

    }

  }
  for (j=1;j<=mesh.kkpll;j++){
    fprintf(out,"%10d",j);
    for (n=1;n<=12;n++)
      fprintf(out,"%10d",vflag[findex(12,j,n)]);
    fprintf(out,"\n");
  }

  delete[] vflag;

  //for each point determine a list of corners
  vflag=new int[8*mesh.kkpll];
  
  for (i=1;i<=mesh.kkpll;i++)
    for (n=1;n<=8;n++)
      vflag[findex(8,i,n)]=-1;
  
  for (i=1;i<=mesh.kkcll;i++) {
    
    z=mesh.kkcz[findex1(i)];

    if (mesh.kztyp[findex1(z)]>0) {

      p1=mesh.kkcp[findex1(i)];
      j=1;
      while (vflag[findex(8,p1,j)]!=-1)
	j=j+1;
      vflag[findex(8,p1,j)]=i;       
    }
  }
  for (j=1;j<=mesh.kkpll;j++){
    fprintf(out,"%10d",j);
    for (n=1;n<=8;n++)
      fprintf(out,"%10d",vflag[findex(8,j,n)]);
    fprintf(out,"\n");
  }
        
  //for each point determine a list of cells
  for (j=1;j<=mesh.kkpll;j++){
    
    for (n=1;n<=8;n++)
      list8[findex1(n)]=-1;
    
    for (i=1;i<=8;i++) {
      k=vflag[findex(8,j,i)];
      if (k!=-1) {
                
	list8[findex1(i)]=mesh.kkcz[findex1(k)];

	if (list8[findex1(i)]>mesh.kkzln) printf("z= %d\n",list8[findex1(i)]);

      }
    }
    fprintf(out,"%10d",j);
    for (n=1;n<=8;n++)
      fprintf(out,"%10d",list8[findex1(n)]);
    fprintf(out,"\n");
  }           

  delete[] vflag;

  //write number of edges
  fprintf(out,"%10d\n",mesh.kkell);

  //for each edge print the 2 points
  for (j=1;j<=mesh.kkell;j++)
    fprintf(out,"%10d%10d%10d\n",j,mesh.kkep1[findex1(j)],mesh.kkep2[findex1(j)]);

  //for each edge determine the list of faces
  vflag=new int[4*mesh.kkell];
  
  for (i=1;i<=mesh.kkell;i++)
    for (n=1;n<=4;n++)
      vflag[findex(4,i,n)]=-1;
  
  for (i=1;i<=mesh.kkfll;i++) {   
    //cell corresponding to this face
    z=mesh.kkfz[findex1(i)];

    if (z>mesh.kkzln) printf("z= %d\n",z);

    //first side in the cell
    j1=mesh.kkzs[findex1(z)]-23;

    //last side in the cell
    j2=mesh.kkzs[findex1(z)];

    //loop over sides in the cell 
    for (j=j1;j<=j2;j++) 
      if (mesh.kksf[findex1(j)]==i) { //this side belongs to the face

	k=mesh.kkse[findex1(j)];
	m=1;
        while (vflag[findex(4,k,m)]!=-1)
	  m=m+1;
 
	//tag edge k with face i
	vflag[findex(4,k,m)]=i;
      }
  }
  for (j=1;j<=mesh.kkell;j++){
    fprintf(out,"%10d",j);
    for (n=1;n<=4;n++)
      fprintf(out,"%10d",vflag[findex(4,j,n)]);
    fprintf(out,"\n");
  }

  delete[] vflag;


  //write number of faces
  fprintf(out,"%10d\n",mesh.kkfll);

  //for each face determine the list of 4 points
  for (i=1;i<=mesh.kkfll;i++) {
 
    //list of points for this face
    for (k=1;k<=4;k++)
      mlist[findex1(k)]=-1;
          
    //cell corresponding to this face
    z=mesh.kkfz[findex1(i)];
 
    //first side in the cell
    j1=mesh.kkzs[findex1(z)]-23;

    //last side in the cell
    j2=mesh.kkzs[findex1(z)];

    //loop over sides in the cell 
    for (j=j1;j<=j2;j++) {
      if (mesh.kksf[findex1(j)]==i) {//this side belongs to the face

	//assign p1 to the face
	p1=mesh.kksp1[findex1(j)];

	idup=0;

	k=1;
	while (mlist[findex1(k)]!=-1) {

	  if (mlist[findex1(k)]==p1) idup=1;
                      
	  k=k+1;

	  if (k==5) break;
	}
                   
	if (idup==0 && k<=4) mlist[findex1(k)]=p1;
                
	//assign p2 to the face
	p2=mesh.kksp2[findex1(j)];

	idup=0;

	k=1;
	while (mlist[findex1(k)]!=-1) {

	  if (mlist[findex1(k)]==p2) idup=1;
                      
	  k=k+1;

	  if (k==5) break;
	}	

	if (idup==0 && k<=4) mlist[findex1(k)]=p2;

      } //end if                     
    }

    fprintf(out,"%10d",i);
    for (n=1;n<=4;n++)
      fprintf(out,"%10d",mlist[findex1(n)]);
    fprintf(out,"\n");
  }
   
  //for each face determine the list of 4 edges
  for (i=1;i<=mesh.kkfll;i++) {
 
    //list of points for this face
    for (k=1;k<=4;k++)
      mlist[findex1(k)]=-1;
          
    //cell corresponding to this face
    z=mesh.kkfz[findex1(i)];
 
    //first side in the cell
    j1=mesh.kkzs[findex1(z)]-23;

    //last side in the cell
    j2=mesh.kkzs[findex1(z)];

    //loop over sides in the cell 
    for (j=j1;j<=j2;j++) {
      if (mesh.kksf[findex1(j)]==i) {//this side belongs to the face
  
	k=mesh.kkse[findex1(j)];
	m=1;
        while (mlist[findex1(m)]!=-1)
	  m=m+1;

	//put edge k in the list
	mlist[findex1(m)]=k;
      }
    }
                
    fprintf(out,"%10d",i);
    for (n=1;n<=4;n++)
      fprintf(out,"%10d",mlist[findex1(n)]);
    fprintf(out,"\n");
  }

  //for each face determine the list of cells on the face
  for (i=1;i<=mesh.kkfll;i++) {
 
    //list of points for this face
    for (k=1;k<=4;k++)
      mlist[findex1(k)]=-1;  
          
    mlist[findex1(1)]=mesh.kkfz[findex1(i)];

    z=mesh.kkfz2[findex1(i)];

    if (z<=mesh.kkzln) mlist[findex1(2)]=z;

    fprintf(out,"%10d%10d%10d\n",i,mlist[findex1(1)],mlist[findex1(2)]);
          
  }


  //write the number of cells
  fprintf(out,"%10d\n",mesh.kkzln);

  //for each cell determine the list of 8 vertices     
  for (i=1;i<=mesh.kkzln;i++) {
    fprintf(out,"%10d",i);
    for (n=1;n<=8;n++)
      fprintf(out,"%10d",mesh.indx[findex(8,i,n)]);
    fprintf(out,"\n");
  }    

 
  //for each cell determine the list of 8 corners
  j=1;
  for (i=1;i<=mesh.kkzln;i++) {
    //make a ordered list of corners in the cell, following order of vertices
    for (k=1;k<=8;k++)
      list8[findex1(k)]=-1;
    
    for (k=1;k<=8;k++)
      for (ii=1;ii<=8;ii++) {  
	m=j+ii-1;
	if (mesh.kkcp[findex1(m)]==mesh.indx[findex(8,i,k)]) list8[findex1(k)]=m;
      }
    
    fprintf(out,"%10d",i);
    for (n=1;n<=8;n++)
      fprintf(out,"%10d",list8[findex1(n)]);
    fprintf(out,"\n");
    
    j=j+8;
  }

  //for each cell determine the list of 6 faces
  for (i=1;i<=mesh.kkzln;i++) {

    for (k=1;k<=8;k++)
      list8[findex1(k)]=-1;  

    //first side in the cell
    j=mesh.kkzs[findex1(i)]-23;

    for (k=1;k<=6;k++) {

      list8[findex1(k)]=mesh.kksf[findex1(j)];

      j=j+4;
    }
    
    fprintf(out,"%10d",i);
    for (n=1;n<=6;n++)
      fprintf(out,"%10d",list8[findex1(n)]);
    fprintf(out,"\n");

  }

  //write the number of corners
  fprintf(out,"%10d\n",8*mesh.kkzln);
       
  //for each corner print the point and cell
  for (i=1;i<=8*mesh.kkzln;i++) 
    fprintf(out,"%10d%10d%10d\n",i, mesh.kkcp[findex1(i)],mesh.kkcz[findex1(i)]);

   
  fclose(out);


}
