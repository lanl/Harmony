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
#include "mpi.h"


class element_t {
  
public :
  int points[8];
  int corners[8];
  int faces[6];
  int flag;
  int cid;

  element_t();
};


element_t::element_t()
{
  int i;

  for (i=0;i<8;i++) {
    points[i]=0;
    corners[i]=0;
  }

  for (i=0;i<6;i++) {  
    faces[i]=0;
  }

  flag=0;
  cid=0;

    
}


class point_t {
  
public :
  double px[3];

  int edges[6];
  int corners[8];
  int elements[8];
  int faces[12];  
  int gid;
  
  point_t();


};

point_t::point_t()
{
  int i;
  
  for (i=0;i<3;i++) 
    px[i]=0.0;

  for (i=0;i<6;i++)
    edges[i]=-1;

  for (i=0;i<8;i++) {
    corners[i]=-1;
    elements[i]=-1;
  }

  for (i=0;i<12;i++)
    faces[i]=-1;

  gid=0;

}


class edge_t {
  
public :
  int points[2];
  int faces[4];  

  edge_t();
};

edge_t::edge_t()
{
  int i;
  
  for (i=0;i<2;i++) 
    points[i]=0;
  
  for (i=0;i<4;i++) 
    faces[i]=-1;  

}

class face_t {
  
public :
  int points[4];
  int edges[4];
  int elements[2];

  face_t();
};

face_t::face_t()
{
  int i;
  
  for (i=0;i<4;i++) { 
    points[i]=0;
    edges[i]=0;
  }
  
  for (i=0;i<2;i++)    
    elements[i]=-1;
}

class corner_t {
  
public :
  int point;
  int element;

  corner_t();
};

corner_t::corner_t()
{
  point=0;
  element=0;
}



class mesh_t {
  
public :
  
  std::vector<edge_t> edges;
  std::vector<corner_t> corners;
  std::vector<face_t> faces;
  std::vector<point_t> points;
  std::vector<element_t> elements;
  

  void make_mesh(char* file);  
  mesh_t(int rank);

  std::vector<int> edges_of_point(int v);
  std::vector<int> corners_of_point(int v);
  std::vector<int> faces_of_point(int v);
  std::vector<int> cells_of_point(int v);

  std::vector<int> vertices();
  std::vector<int> cells_owned();
  std::vector<int> cells();
  std::vector<int> vertices_of_cell(int c);
  std::vector<int> corners_of_cell(int c);
  std::vector<int> faces_of_cell(int c);
  std::vector<int> list_faces();
  std::vector<int> list_edges();  
  std::vector<int> list_corners();

  std::vector<int> edges_of_face(int f);
  std::vector<int> points_of_face(int f);
  std::vector<int> cells_of_face(int f);
  std::vector<int> points_of_edge(int e);
  std::vector<int> faces_of_edge(int e);
  
};

std::vector<int> mesh_t::edges_of_point(int v)
{

  std::vector<int> result;

  result.reserve(6);

  auto &p = points.at(v);
  
  for (int i=0;i<6;i++)
    if (p.edges[i]>-1) result.push_back(p.edges[i]);

  return result;
  
}

std::vector<int> mesh_t::corners_of_point(int v)
{

  std::vector<int> result;

  result.reserve(8);
  
  auto &p = points.at(v);
  
  for (int i=0;i<8;i++)
    if (p.corners[i]>-1) result.push_back(p.corners[i]);

  return result;
  
}

std::vector<int> mesh_t::cells_of_point(int v)
{

  std::vector<int> result;

  result.reserve(8);
  
  auto &p = points.at(v);
  
  for (int i=0;i<8;i++)
    if (p.elements[i]>-1) result.push_back(p.elements[i]);

  return result;
  
}

std::vector<int> mesh_t::faces_of_point(int v)
{

  std::vector<int> result;

  result.reserve(12);
  
  auto &p = points.at(v);
  
  for (int i=0;i<12;i++)
    if (p.faces[i]>-1) result.push_back(p.faces[i]);

  return result;
  
}

std::vector<int> mesh_t::vertices()
{

  std::vector<int> result;

  int num=points.size();

  result.reserve(num);

  for (int i=0;i<num;i++) 
    result.push_back(i);

  return result;
  
}

std::vector<int> mesh_t::cells_owned()
{

  std::vector<int> result;

  int num=elements.size();

  result.reserve(num);

  for (int i=0;i<num;i++) {
    
    auto &c = elements.at(i);
    if (c.flag) result.push_back(i);

  }

  return result;
  
}


std::vector<int> mesh_t::cells()
{

  std::vector<int> result;

  int num=elements.size();

  result.reserve(num);

  for (int i=0;i<num;i++) 
    result.push_back(i);

  return result;
  
}

std::vector<int> mesh_t::vertices_of_cell(int c)
{

  std::vector<int> result;

  result.reserve(8);
    
  auto &ic = elements.at(c);

  for (int i=0;i<8;i++)
    result.push_back(ic.points[i]);

  return result;
  
}

std::vector<int> mesh_t::corners_of_cell(int c)
{

  std::vector<int> result;

  result.reserve(8);
    
  auto &ic = elements.at(c);

  for (int i=0;i<8;i++)
    result.push_back(ic.corners[i]);

  return result;
  
}

std::vector<int> mesh_t::faces_of_cell(int c)
{

  std::vector<int> result;

  result.reserve(6);
    
  auto &ic = elements.at(c);

  for (int i=0;i<6;i++)
    result.push_back(ic.faces[i]);

  return result;
  
}



std::vector<int> mesh_t::edges_of_face(int f)
{

  std::vector<int> result;

  result.reserve(4);
    
  auto &fc = faces.at(f);

  for (int i=0;i<4;i++)
    result.push_back(fc.edges[i]);

  return result;
  
}

std::vector<int> mesh_t::points_of_face(int f)
{

  std::vector<int> result;

  result.reserve(4);
    
  auto &fc = faces.at(f);

  for (int i=0;i<4;i++)
    result.push_back(fc.points[i]);

  return result;
  
}

std::vector<int> mesh_t::cells_of_face(int f)
{

  std::vector<int> result;

  result.reserve(2);
    
  auto &fc = faces.at(f);

  for (int i=0;i<2;i++)
    if (fc.elements[i]>-1) result.push_back(fc.elements[i]);

  return result;
  
}

std::vector<int> mesh_t::points_of_edge(int e)
{

  std::vector<int> result;

  result.reserve(2);
    
  auto &ec = edges.at(e);

  for (int i=0;i<2;i++)
    result.push_back(ec.points[i]);

  return result;
  
}

std::vector<int> mesh_t::list_faces()
{

  std::vector<int> result;

  int num=faces.size();

  result.reserve(num);

  for (int i=0;i<num;i++) 
    result.push_back(i);

  return result;
  
}

std::vector<int> mesh_t::list_edges()
{

  std::vector<int> result;

  int num=edges.size();

  result.reserve(num);

  for (int i=0;i<num;i++) 
    result.push_back(i);

  return result;
  
}

std::vector<int> mesh_t::list_corners()
{

  std::vector<int> result;

  int num=corners.size();

  result.reserve(num);

  for (int i=0;i<num;i++) 
    result.push_back(i);

  return result;
  
}

std::vector<int> mesh_t::faces_of_edge(int e)
{

  std::vector<int> result;

  result.reserve(4);
    
  auto &ec = edges.at(e);

  for (int i=0;i<4;i++)
    if (ec.faces[i]>-1) result.push_back(ec.faces[i]);

  return result;
  
}

void mesh_t::make_mesh(char* file)  
{

  int i,j,num,indx;

  char c;
  
  FILE *in;

  int kkpll,kkzln;

  double xa;

  in=fopen(file,"r");

  //read vertex information
  fscanf(in,"%d",&kkpll);

  points.reserve(kkpll);
  points.resize(kkpll);

  printf("size of vertex array= %d\n",points.size());  

  for (i=0;i<kkpll;i++) {

    //read global id for the node
    fscanf(in,"%d",&num);
    points[i].gid=num;

    //read vertex positions
    fscanf(in,"%le",&xa);
    points.at(i).px[0]=xa;

    fscanf(in,"%le",&xa);
    points[i].px[1]=xa;

    fscanf(in,"%le",&xa);
    points[i].px[2]=xa;

    //read edge connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      points[i].edges[j]=indx;
    }

    //read corner connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      points[i].corners[j]=indx;
    }
    
    //read face connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      points[i].faces[j]=indx;
    }

    //read element connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      points[i].elements[j]=indx;
    }    

    
  } //end vertex setup


  //read edge information
  fscanf(in,"%d",&kkpll);
  edges.reserve(kkpll);
  edges.resize(kkpll);

  printf("size of edge array= %d\n",edges.size());  

  for (i=0;i<kkpll;i++) {

    //read point connectivities
    fscanf(in,"%d",&num);
     
    for (j=0;j<2;j++) {
      fscanf(in,"%d",&indx);
      edges.at(i).points[j]=indx;
    }
    
    //read face connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      edges[i].faces[j]=indx;
    }

  } //done reading edge information


  //read face information
  fscanf(in,"%d",&kkpll);
  faces.reserve(kkpll);
  faces.resize(kkpll);

  printf("size of face array= %d\n",faces.size());  

  for (i=0;i<kkpll;i++) {

    //read point connectivities
    fscanf(in,"%d",&num);
    
    for (j=0;j<4;j++) {
      fscanf(in,"%d",&indx);
      faces.at(i).points[j]=indx;
    }
    
    //read edge connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      faces[i].edges[j]=indx;
    }

    //read element connectivities
    fscanf(in,"%d",&num);

    for (j=0;j<num;j++) {
      fscanf(in,"%d",&indx);
      faces[i].elements[j]=indx;
    }    

  } //done reading face information  
    

  //read element information
  fscanf(in,"%d",&kkpll);
  elements.reserve(kkpll);
  elements.resize(kkpll);

  printf("size of element array= %d\n",elements.size());  

  for (i=0;i<kkpll;i++) {

    //read global id for the cell
    fscanf(in,"%d",&num);
    elements[i].cid=num;

    //read ownership flag
    fscanf(in,"%d",&num);
    elements.at(i).flag=num;

    //read point connectivities
    fscanf(in,"%d",&num);   
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      elements[i].points[j]=indx;
    }
    
    //read corner connectivities
    fscanf(in,"%d",&num);       
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      elements[i].corners[j]=indx;
    }

    //read face connectivities
    fscanf(in,"%d",&num);       
    for (j=0;j<6;j++) {
      fscanf(in,"%d",&indx);
      elements[i].faces[j]=indx;
    }    

  } //done reading element information
  
  //read corner information
  fscanf(in,"%d",&kkpll);
  corners.reserve(kkpll);
  corners.resize(kkpll);

  printf("size of corner array= %d\n",corners.size());  

  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&indx);
    corners[i].point=indx;
    
    fscanf(in,"%d",&indx);
    corners[i].element=indx; 

    // printf("corner = %d node = %d cell = %d\n",i, corners[i].point,corners[i].element);  
  } //done reading corner information

  fclose(in);

}

mesh_t::mesh_t(int rank)  
{

  int i,j,num,indx;

  char c;
  
  FILE *in;

  int kkpll,kkzln;

  double xa;

  char file_name[100];

  char ch;
  
  sprintf(file_name,"mesh/ensight_%05d/mesh_file.dat",rank);  

  in=fopen(file_name,"r");

  //read vertex information
  fscanf(in,"%d",&kkpll);
  points.reserve(kkpll);
  points.resize(kkpll);

  printf("size of vertex array= %d\n",points.size());
  
 
  //read edge connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<6;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      points[i].edges[j]=indx;
    }
  }
    
  //read face connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<12;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      points[i].faces[j]=indx;
    }
  }

  //read corner connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      points[i].corners[j]=indx;
    }

  }  

  //read element connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    // if (rank==0) printf("point = %d ",i);
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      points[i].elements[j]=indx;
      // if (rank==0) printf("%d ",points[i].elements[j]);
    }
    //if (rank==0) printf("\n");
  } //end vertex setup


  //read edge information
  fscanf(in,"%d",&kkpll);
  edges.reserve(kkpll);
  edges.resize(kkpll);

  printf("size of edge array= %d\n",edges.size());
  
  //read point connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<2;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      edges[i].points[j]=indx;
    }
  }
    
  //read face connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<4;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      edges[i].faces[j]=indx;
    }

  } //done reading edge information


  //read face information
  fscanf(in,"%d",&kkpll);
  faces.reserve(kkpll);
  faces.resize(kkpll);

  printf("size of face array= %d\n",faces.size());  

  //read point connectivities  
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<4;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      faces[i].points[j]=indx;
    }
  }
    
  //read edge connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<4;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      faces[i].edges[j]=indx;
    }
  }

  //read element connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<2;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      faces[i].elements[j]=indx;
    }    
  } //done reading face information  
    

  //read element information
  fscanf(in,"%d",&kkpll);
  elements.reserve(kkpll);
  elements.resize(kkpll);

  printf("size of element array= %d\n",elements.size());

  
  //read point connectivities
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      elements[i].points[j]=indx;
    }
  }
    
 
 //read corner connectivities
 for (i=0;i<kkpll;i++) {
   fscanf(in,"%d",&num);
    for (j=0;j<8;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      elements[i].corners[j]=indx;
    }
 }

 //read face connectivities
 for (i=0;i<kkpll;i++) {
   fscanf(in,"%d",&num);
    for (j=0;j<6;j++) {
      fscanf(in,"%d",&indx);
      if (indx!=-1) indx=indx-1;
      elements[i].faces[j]=indx;
    }    

  }

 //done reading element information
  
  //read corner information
  fscanf(in,"%d",&kkpll);
  corners.reserve(kkpll);
  corners.resize(kkpll);

  printf("size of corner array= %d\n",corners.size());  

  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    
    fscanf(in,"%d",&indx);
    corners[i].point=indx-1;
    
    fscanf(in,"%d",&indx);
    corners[i].element=indx-1; 

    // printf("corner = %d node = %d cell = %d\n",i, corners[i].point,corners[i].element);  
  } //done reading corner information

  fclose(in);


  sprintf(file_name,"mesh/ensight_%05d/engold.geo",rank);
  printf("%s\n",file_name);
  
  in=fopen(file_name,"r");

  for (j=1;j<=12;j++) {

    i=0;
    while ((ch=(char)fgetc(in))!='\n') {
        i++;
    }

  }

  fscanf(in,"%d",&j);
  printf("%d\n",j);

  kkpll=j;
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&j);
  }
  
  
  //read vertex positions
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%le",&xa);
    points[i].px[0]=xa;
  }
  
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%le",&xa);
    points[i].px[1]=xa;
  }
  
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%le",&xa);
    points[i].px[2]=xa;
  }
  
  ch=(char)fgetc(in);
  for (j=1;j<=1;j++) {

    char buffer[100];

    i=0;
    while ((ch=(char)fgetc(in))!='\n') {
      buffer[i]=ch;
        i++;
    }
    buffer[i]='\0';
    printf("rank= %d buf= %s\n",rank,buffer);

  }

  fscanf(in,"%d",&j);
  printf("%d\n",j);

  kkpll=j;
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&j);
  }

  int pn[8];
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d%d%d%d%d%d%d%d",pn,pn+1,pn+2,pn+3,pn+4,pn+5,pn+6,pn+7);

    for (j=0;j<8;j++) {
      pn[j]--;
      if (pn[j]!=elements[i].points[j]) printf("rank= %d cell=%d point order not correct\n",rank,i);
    }
    
  }  
  

  fclose(in);

  sprintf(file_name,"mesh/ensight_%05d/engold.flg",rank);  
  printf("%s\n",file_name);
  
  in=fopen(file_name,"r");

  kkpll=elements.size();
  
  //read cell ownership flag
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    elements[i].flag=num;
  }

  fclose(in);

  sprintf(file_name,"mesh/ensight_%05d/engold.cid",rank);  
  printf("%s\n",file_name);
  
  in=fopen(file_name,"r");

  kkpll=elements.size();
  
  //read cell global id
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    elements[i].cid=num;
  }

  fclose(in);  

  sprintf(file_name,"mesh/ensight_%05d/engold.gid",rank);  
  printf("%s\n",file_name);
  
  in=fopen(file_name,"r");

  kkpll=points.size();
  
  //read point global id
  for (i=0;i<kkpll;i++) {
    fscanf(in,"%d",&num);
    points[i].gid=num;
  }

  fclose(in);    
  
  
}



