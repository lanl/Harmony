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
typedef struct _mydata2 { /* data type used for computing gradients in "slopes2" and "block2" */

  double vol,sx[3],sy[3],sz[3];
  int p[8];
  double u[27],v[27],w[27];

} mydata2_t; 
 

void block2(mymesh_t & mesh,mydata2_t & mydata) /*subroutine used by "slopes2" to find node-centered gradients*/
{
  int i,j,k;
  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  double b11,b12,b13,b14,b15,b16,b17,b18;
  double b21,b22,b23,b24,b25,b26,b27,b28;
  double b31,b32,b33,b34,b35,b36,b37,b38;

  /*centroids of the neighboring vertices forming the block 2 control volume*/

  x1=mesh.xc[mydata.p[0]];
  x2=mesh.xc[mydata.p[1]];
  x3=mesh.xc[mydata.p[2]];
  x4=mesh.xc[mydata.p[3]];
  x5=mesh.xc[mydata.p[4]];
  x6=mesh.xc[mydata.p[5]];
  x7=mesh.xc[mydata.p[6]];
  x8=mesh.xc[mydata.p[7]];

  y1=mesh.yc[mydata.p[0]];
  y2=mesh.yc[mydata.p[1]];
  y3=mesh.yc[mydata.p[2]];
  y4=mesh.yc[mydata.p[3]];
  y5=mesh.yc[mydata.p[4]];
  y6=mesh.yc[mydata.p[5]];
  y7=mesh.yc[mydata.p[6]];
  y8=mesh.yc[mydata.p[7]];

  z1=mesh.zc[mydata.p[0]];
  z2=mesh.zc[mydata.p[1]];
  z3=mesh.zc[mydata.p[2]];
  z4=mesh.zc[mydata.p[3]];
  z5=mesh.zc[mydata.p[4]];
  z6=mesh.zc[mydata.p[5]];
  z7=mesh.zc[mydata.p[6]];
  z8=mesh.zc[mydata.p[7]];


  /* volume of block 2 */

  mydata.vol+=(x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
                   x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
                   x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
                   x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
                   x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
                   x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
                   x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
	       x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;


  /* B matrix */

  b11=(-y2*z3 -y2*z4 +y2*z5 +y2*z6 +y3*z2 -y3*z4 +y4*z2 +y4*z3 -y4*z5 -y4*z8 -y5*z2 +y5*z4 -y5*z6 +y5*z8 -y6*z2 +y6*z5 +y8*z4 -y8*z5 )/12.;

  b12=(+y1*z3 +y1*z4 -y1*z5 -y1*z6 -y3*z1 -y3*z4 +y3*z6 +y3*z7 -y4*z1 +y4*z3 +y5*z1 -y5*z6 +y6*z1 -y6*z3 +y6*z5 -y6*z7 -y7*z3 +y7*z6 )/12.;

  b13=(-y1*z2 +y1*z4 +y2*z1 +y2*z4 -y2*z6 -y2*z7 -y4*z1 -y4*z2 +y4*z7 +y4*z8 +y6*z2 -y6*z7 +y7*z2 -y7*z4 +y7*z6 -y7*z8 -y8*z4 +y8*z7 )/12.;

  b14=(-y1*z2 -y1*z3 +y1*z5 +y1*z8 +y2*z1 -y2*z3 +y3*z1 +y3*z2 -y3*z7 -y3*z8 -y5*z1 +y5*z8 +y7*z3 -y7*z8 -y8*z1 +y8*z3 -y8*z5 +y8*z7 )/12.;

  b15=(+y1*z2 -y1*z4 +y1*z6 -y1*z8 -y2*z1 +y2*z6 +y4*z1 -y4*z8 -y6*z1 -y6*z2 +y6*z7 +y6*z8 -y7*z6 +y7*z8 +y8*z1 +y8*z4 -y8*z6 -y8*z7 )/12.;

  b16=(+y1*z2 -y1*z5 -y2*z1 +y2*z3 -y2*z5 +y2*z7 -y3*z2 +y3*z7 +y5*z1 +y5*z2 -y5*z7 -y5*z8 -y7*z2 -y7*z3 +y7*z5 +y7*z8 +y8*z5 -y8*z7)/12.;

  b17=(+y2*z3 -y2*z6 -y3*z2 +y3*z4 -y3*z6 +y3*z8 -y4*z3 +y4*z8 +y5*z6 -y5*z8 +y6*z2 +y6*z3 -y6*z5 -y6*z8 -y8*z3 -y8*z4 +y8*z5 +y8*z6 )/12.;

  b18=(-y1*z4 +y1*z5 +y3*z4 -y3*z7 +y4*z1 -y4*z3 +y4*z5 -y4*z7 -y5*z1 -y5*z4 +y5*z6 +y5*z7 -y6*z5 +y6*z7 +y7*z3 +y7*z4 -y7*z5 -y7*z6 )/12.;


  b21=(-z2*x3 -z2*x4 +z2*x5 +z2*x6 +z3*x2 -z3*x4 +z4*x2 +z4*x3 -z4*x5 -z4*x8 -z5*x2 +z5*x4 -z5*x6 +z5*x8 -z6*x2 +z6*x5 +z8*x4 -z8*x5 )/12.;

  b22=(+z1*x3 +z1*x4 -z1*x5 -z1*x6 -z3*x1 -z3*x4 +z3*x6 +z3*x7 -z4*x1 +z4*x3 +z5*x1 -z5*x6 +z6*x1 -z6*x3 +z6*x5 -z6*x7 -z7*x3 +z7*x6 )/12.;

  b23=(-z1*x2 +z1*x4 +z2*x1 +z2*x4 -z2*x6 -z2*x7 -z4*x1 -z4*x2 +z4*x7 +z4*x8 +z6*x2 -z6*x7 +z7*x2 -z7*x4 +z7*x6 -z7*x8 -z8*x4 +z8*x7 )/12.;

  b24=(-z1*x2 -z1*x3 +z1*x5 +z1*x8 +z2*x1 -z2*x3 +z3*x1 +z3*x2 -z3*x7 -z3*x8 -z5*x1 +z5*x8 +z7*x3 -z7*x8 -z8*x1 +z8*x3 -z8*x5 +z8*x7 )/12.;

  b25=(+z1*x2 -z1*x4 +z1*x6 -z1*x8 -z2*x1 +z2*x6 +z4*x1 -z4*x8 -z6*x1 -z6*x2 +z6*x7 +z6*x8 -z7*x6 +z7*x8 +z8*x1 +z8*x4 -z8*x6 -z8*x7 )/12.;

  b26=(+z1*x2 -z1*x5 -z2*x1 +z2*x3 -z2*x5 +z2*x7 -z3*x2 +z3*x7 +z5*x1 +z5*x2 -z5*x7 -z5*x8 -z7*x2 -z7*x3 +z7*x5 +z7*x8 +z8*x5 -z8*x7 )/12.;

  b27=(+z2*x3 -z2*x6 -z3*x2 +z3*x4 -z3*x6 +z3*x8 -z4*x3 +z4*x8 +z5*x6 -z5*x8 +z6*x2 +z6*x3 -z6*x5 -z6*x8 -z8*x3 -z8*x4 +z8*x5 +z8*x6 )/12.;

  b28=(-z1*x4 +z1*x5 +z3*x4 -z3*x7 +z4*x1 -z4*x3 +z4*x5 -z4*x7 -z5*x1 -z5*x4 +z5*x6 +z5*x7 -z6*x5 +z6*x7 +z7*x3 +z7*x4 -z7*x5 -z7*x6 )/12.;


  b31=(-x2*y3 -x2*y4 +x2*y5 +x2*y6 +x3*y2 -x3*y4 +x4*y2 +x4*y3 -x4*y5 -x4*y8 -x5*y2 +x5*y4 -x5*y6 +x5*y8 -x6*y2 +x6*y5 +x8*y4 -x8*y5 )/12.;

  b32=(+x1*y3 +x1*y4 -x1*y5 -x1*y6 -x3*y1 -x3*y4 +x3*y6 +x3*y7 -x4*y1 +x4*y3 +x5*y1 -x5*y6 +x6*y1 -x6*y3 +x6*y5 -x6*y7 -x7*y3 +x7*y6 )/12.;

  b33=(-x1*y2 +x1*y4 +x2*y1 +x2*y4 -x2*y6 -x2*y7 -x4*y1 -x4*y2 +x4*y7 +x4*y8 +x6*y2 -x6*y7 +x7*y2 -x7*y4 +x7*y6 -x7*y8 -x8*y4 +x8*y7 )/12.;

  b34=(-x1*y2 -x1*y3 +x1*y5 +x1*y8 +x2*y1 -x2*y3 +x3*y1 +x3*y2 -x3*y7 -x3*y8 -x5*y1 +x5*y8 +x7*y3 -x7*y8 -x8*y1 +x8*y3 -x8*y5 +x8*y7 )/12.;

  b35=(+x1*y2 -x1*y4 +x1*y6 -x1*y8 -x2*y1 +x2*y6 +x4*y1 -x4*y8 -x6*y1 -x6*y2 +x6*y7 +x6*y8 -x7*y6 +x7*y8 +x8*y1 +x8*y4 -x8*y6 -x8*y7 )/12.;

  b36=(+x1*y2 -x1*y5 -x2*y1 +x2*y3 -x2*y5 +x2*y7 -x3*y2 +x3*y7 +x5*y1 +x5*y2 -x5*y7 -x5*y8 -x7*y2 -x7*y3 +x7*y5 +x7*y8 +x8*y5 -x8*y7 )/12.;

  b37=(+x2*y3 -x2*y6 -x3*y2 +x3*y4 -x3*y6 +x3*y8 -x4*y3 +x4*y8 +x5*y6 -x5*y8 +x6*y2 +x6*y3 -x6*y5 -x6*y8 -x8*y3 -x8*y4 +x8*y5 +x8*y6 )/12.;

  b38=(-x1*y4 +x1*y5 +x3*y4 -x3*y7 +x4*y1 -x4*y3 +x4*y5 -x4*y7 -x5*y1 -x5*y4 +x5*y6 +x5*y7 -x6*y5 +x6*y7 +x7*y3 +x7*y4 -x7*y5 -x7*y6 )/12.;


  /* block 2 u velocity derivative */
  mydata.sx[0]+=(mydata.u[mydata.p[0]]*b11+mydata.u[mydata.p[1]]*b12+mydata.u[mydata.p[2]]*b13+mydata.u[mydata.p[3]]*b14+mydata.u[mydata.p[4]]*b15+mydata.u[mydata.p[5]]*b16+mydata.u[mydata.p[6]]*b17+mydata.u[mydata.p[7]]*b18);
  mydata.sy[0]+=(mydata.u[mydata.p[0]]*b21+mydata.u[mydata.p[1]]*b22+mydata.u[mydata.p[2]]*b23+mydata.u[mydata.p[3]]*b24+mydata.u[mydata.p[4]]*b25+mydata.u[mydata.p[5]]*b26+mydata.u[mydata.p[6]]*b27+mydata.u[mydata.p[7]]*b28);
  mydata.sz[0]+=(mydata.u[mydata.p[0]]*b31+mydata.u[mydata.p[1]]*b32+mydata.u[mydata.p[2]]*b33+mydata.u[mydata.p[3]]*b34+mydata.u[mydata.p[4]]*b35+mydata.u[mydata.p[5]]*b36+mydata.u[mydata.p[6]]*b37+mydata.u[mydata.p[7]]*b38);

  /* block 2 v velocity derivative */
  mydata.sx[1]+=(mydata.v[mydata.p[0]]*b11+mydata.v[mydata.p[1]]*b12+mydata.v[mydata.p[2]]*b13+mydata.v[mydata.p[3]]*b14+mydata.v[mydata.p[4]]*b15+mydata.v[mydata.p[5]]*b16+mydata.v[mydata.p[6]]*b17+mydata.v[mydata.p[7]]*b18);
  mydata.sy[1]+=(mydata.v[mydata.p[0]]*b21+mydata.v[mydata.p[1]]*b22+mydata.v[mydata.p[2]]*b23+mydata.v[mydata.p[3]]*b24+mydata.v[mydata.p[4]]*b25+mydata.v[mydata.p[5]]*b26+mydata.v[mydata.p[6]]*b27+mydata.v[mydata.p[7]]*b28);
  mydata.sz[1]+=(mydata.v[mydata.p[0]]*b31+mydata.v[mydata.p[1]]*b32+mydata.v[mydata.p[2]]*b33+mydata.v[mydata.p[3]]*b34+mydata.v[mydata.p[4]]*b35+mydata.v[mydata.p[5]]*b36+mydata.v[mydata.p[6]]*b37+mydata.v[mydata.p[7]]*b38);

  /* block 2 w velocity derivative */
  mydata.sx[2]+=(mydata.w[mydata.p[0]]*b11+mydata.w[mydata.p[1]]*b12+mydata.w[mydata.p[2]]*b13+mydata.w[mydata.p[3]]*b14+mydata.w[mydata.p[4]]*b15+mydata.w[mydata.p[5]]*b16+mydata.w[mydata.p[6]]*b17+mydata.w[mydata.p[7]]*b18);
  mydata.sy[2]+=(mydata.w[mydata.p[0]]*b21+mydata.w[mydata.p[1]]*b22+mydata.w[mydata.p[2]]*b23+mydata.w[mydata.p[3]]*b24+mydata.w[mydata.p[4]]*b25+mydata.w[mydata.p[5]]*b26+mydata.w[mydata.p[6]]*b27+mydata.w[mydata.p[7]]*b28);
  mydata.sz[2]+=(mydata.w[mydata.p[0]]*b31+mydata.w[mydata.p[1]]*b32+mydata.w[mydata.p[2]]*b33+mydata.w[mydata.p[3]]*b34+mydata.w[mydata.p[4]]*b35+mydata.w[mydata.p[5]]*b36+mydata.w[mydata.p[6]]*b37+mydata.w[mydata.p[7]]*b38);

}

void slopes2( 
  mesh_t  &mesh,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {
  /* finds vertex centered derivates of u,v,w velocities */

  double factor,phi[4],fmin[4],fmax[4],xc[8],yc[8],zc[8];

  int m,n,l,i,j,k,kk,ii,jj;

  int flag;

  mydata2_t mydata;
  mymesh_t mymesh;

  int iorder[8],cn[8],myvert[8];


  //cells
  auto cs = mesh.cells();
  auto verts = mesh.vertices();
  auto NP = verts.size();
  auto cs_owned = mesh.cells_owned();
  
  auto c1 = verts[0];
  auto c2 = verts[0];
  auto c3 = verts[0];
  auto c4 = verts[0];

  int rank,numtasks,kk2,k2;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  
  double *recv,*send;

  MPI_Status stat;

  int *vflag;
  vflag=new int[NP];

  // loop over all vertices and tag all non-ghost vertices
  for (auto v : verts )  {
    i=v;
    vflag[i]=0;
    
  }
  
  for (auto c : cs_owned )  {

    for (auto v : mesh.vertices_of_cell(c) )  {
        i=v;
	vflag[i]=1;
    }  
    
  }      
  
  /* loop over all verticies */
  for ( auto v : verts ) {

   if (vflag[v]==1) { //skip ghost vertices

    //begin construction of the indx array containing vertex neighbors
    i=v;
    for (j=0;j<27;j++)
      mymesh.indx[j]=-1; //default setting of -1 indicates boundary vertex

    mymesh.indx[13]=v; //central vertex on the 27 point stencil

    //loop over all cells containing vertex v
    for ( auto c : mesh.cells_of_point(v) ) {

      //order the vertices in the cell c
      myvert[0]=cell[c].iord1;
      myvert[1]=cell[c].iord2;
      myvert[2]=cell[c].iord3;
      myvert[3]=cell[c].iord4;
      myvert[4]=cell[c].iord5;
      myvert[5]=cell[c].iord6;
      myvert[6]=cell[c].iord7;
      myvert[7]=cell[c].iord8;

      j=0;
      for ( auto v1 : mesh.vertices_of_cell(c) ) {
	for (k=0;k<8;k++) 
	  if (node[v1].gid==myvert[k]) iorder[j]=k;

	j++;
      }

    
      k=0;
      for ( auto v1 : mesh.vertices_of_cell(c) ) {
	cn[iorder[k]]=v1;
	k++;
      }

      //find which vertex corresponds to v
      for (k=0,l=0;k<8;k++) {
        if (cn[k]==i) l=k;
      }

      if (l==0) { 

	mymesh.indx[13]=cn[0];
	mymesh.indx[14]=cn[1];
	mymesh.indx[17]=cn[2];
	mymesh.indx[16]=cn[3];
	mymesh.indx[22]=cn[4];
	mymesh.indx[23]=cn[5];
	mymesh.indx[26]=cn[6];
	mymesh.indx[25]=cn[7];

	mymesh.xc[13]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[14]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[17]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[16]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[22]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[23]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[26]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[25]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[13]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[14]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[17]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[16]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[22]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[23]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[26]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[25]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[13]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[14]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[17]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[16]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[22]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[23]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[26]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[25]=mesh.points[verts[cn[7]]].px[2];
      }
      else if (l==1) { 

	mymesh.indx[12]=cn[0];
	mymesh.indx[13]=cn[1];
	mymesh.indx[16]=cn[2];
	mymesh.indx[15]=cn[3];
	mymesh.indx[21]=cn[4];
	mymesh.indx[22]=cn[5];
	mymesh.indx[25]=cn[6];
	mymesh.indx[24]=cn[7];

	mymesh.xc[12]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[16]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[15]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[21]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[22]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[25]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[24]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[12]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[16]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[15]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[21]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[22]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[25]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[24]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[12]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[16]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[15]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[21]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[22]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[25]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[24]=mesh.points[verts[cn[7]]].px[2];
      }

      else if (l==2) { 

	mymesh.indx[9]=cn[0];
	mymesh.indx[10]=cn[1];
	mymesh.indx[13]=cn[2];
	mymesh.indx[12]=cn[3];
	mymesh.indx[18]=cn[4];
	mymesh.indx[19]=cn[5];
	mymesh.indx[22]=cn[6];
	mymesh.indx[21]=cn[7];

	mymesh.xc[9]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[10]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[12]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[18]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[19]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[22]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[21]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[9]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[10]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[12]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[18]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[19]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[22]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[21]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[9]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[10]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[12]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[18]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[19]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[22]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[21]=mesh.points[verts[cn[7]]].px[2];
      }
      else if (l==3) { 

	mymesh.indx[10]=cn[0];
	mymesh.indx[11]=cn[1];
	mymesh.indx[14]=cn[2];
	mymesh.indx[13]=cn[3];
	mymesh.indx[19]=cn[4];
	mymesh.indx[20]=cn[5];
	mymesh.indx[23]=cn[6];
	mymesh.indx[22]=cn[7];

	mymesh.xc[10]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[11]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[14]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[19]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[20]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[23]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[22]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[10]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[11]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[14]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[19]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[20]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[23]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[22]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[10]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[11]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[14]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[19]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[20]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[23]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[22]=mesh.points[verts[cn[7]]].px[2];
      }

      else if (l==4) { 

	mymesh.indx[4]=cn[0];
	mymesh.indx[5]=cn[1];
	mymesh.indx[8]=cn[2];
	mymesh.indx[7]=cn[3];
	mymesh.indx[13]=cn[4];
	mymesh.indx[14]=cn[5];
	mymesh.indx[17]=cn[6];
	mymesh.indx[16]=cn[7];

	mymesh.xc[4]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[5]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[8]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[7]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[14]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[17]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[16]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[4]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[5]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[8]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[7]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[14]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[17]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[16]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[4]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[5]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[8]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[7]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[14]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[17]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[16]=mesh.points[verts[cn[7]]].px[2];
      }
      else if (l==5) { 

	mymesh.indx[3]=cn[0];
	mymesh.indx[4]=cn[1];
	mymesh.indx[7]=cn[2];
	mymesh.indx[6]=cn[3];
	mymesh.indx[12]=cn[4];
	mymesh.indx[13]=cn[5];
	mymesh.indx[16]=cn[6];
	mymesh.indx[15]=cn[7];

	mymesh.xc[3]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[4]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[7]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[6]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[12]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[16]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[15]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[3]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[4]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[7]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[6]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[12]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[16]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[15]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[3]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[4]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[7]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[6]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[12]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[16]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[15]=mesh.points[verts[cn[7]]].px[2];
      }
      else if (l==6) { 

	mymesh.indx[0]=cn[0];
	mymesh.indx[1]=cn[1];
	mymesh.indx[4]=cn[2];
	mymesh.indx[3]=cn[3];
	mymesh.indx[9]=cn[4];
	mymesh.indx[10]=cn[5];
	mymesh.indx[13]=cn[6];
	mymesh.indx[12]=cn[7];

	mymesh.xc[0]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[1]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[4]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[3]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[9]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[10]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[12]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[0]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[1]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[4]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[3]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[9]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[10]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[12]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[0]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[1]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[4]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[3]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[9]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[10]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[12]=mesh.points[verts[cn[7]]].px[2];
      }
      else if (l==7) { 

	mymesh.indx[1]=cn[0];
	mymesh.indx[2]=cn[1];
	mymesh.indx[5]=cn[2];
	mymesh.indx[4]=cn[3];
	mymesh.indx[10]=cn[4];
	mymesh.indx[11]=cn[5];
	mymesh.indx[14]=cn[6];
	mymesh.indx[13]=cn[7];

	mymesh.xc[1]=mesh.points[verts[cn[0]]].px[0];
	mymesh.xc[2]=mesh.points[verts[cn[1]]].px[0];
	mymesh.xc[5]=mesh.points[verts[cn[2]]].px[0];
	mymesh.xc[4]=mesh.points[verts[cn[3]]].px[0];
	mymesh.xc[10]=mesh.points[verts[cn[4]]].px[0];
	mymesh.xc[11]=mesh.points[verts[cn[5]]].px[0];
	mymesh.xc[14]=mesh.points[verts[cn[6]]].px[0];
	mymesh.xc[13]=mesh.points[verts[cn[7]]].px[0];

	mymesh.yc[1]=mesh.points[verts[cn[0]]].px[1];
	mymesh.yc[2]=mesh.points[verts[cn[1]]].px[1];
	mymesh.yc[5]=mesh.points[verts[cn[2]]].px[1];
	mymesh.yc[4]=mesh.points[verts[cn[3]]].px[1];
	mymesh.yc[10]=mesh.points[verts[cn[4]]].px[1];
	mymesh.yc[11]=mesh.points[verts[cn[5]]].px[1];
	mymesh.yc[14]=mesh.points[verts[cn[6]]].px[1];
	mymesh.yc[13]=mesh.points[verts[cn[7]]].px[1];

	mymesh.zc[1]=mesh.points[verts[cn[0]]].px[2];
	mymesh.zc[2]=mesh.points[verts[cn[1]]].px[2];
	mymesh.zc[5]=mesh.points[verts[cn[2]]].px[2];
	mymesh.zc[4]=mesh.points[verts[cn[3]]].px[2];
	mymesh.zc[10]=mesh.points[verts[cn[4]]].px[2];
	mymesh.zc[11]=mesh.points[verts[cn[5]]].px[2];
	mymesh.zc[14]=mesh.points[verts[cn[6]]].px[2];
	mymesh.zc[13]=mesh.points[verts[cn[7]]].px[2];
      }

    } //finished looping over all cells containing vertex v
    // now all entries in the indx array that are not -1 have been set

    //fill the ghost and boundary vertices on the stencil
    complete_stencil(mymesh);

    //setup mydata
    mydata.vol=0.0;

    for (n=0;n<3;n++) {
      mydata.sx[n]=0.0;
      mydata.sy[n]=0.0;
      mydata.sz[n]=0.0;
    }
    
    //copy velocity onto local mesh
    for (n=0;n<27;n++) {
      
      c1=verts[mymesh.cpy[n]];
      mydata.u[n]=node[c1].un[0];
      mydata.v[n]=node[c1].un[1];
      mydata.w[n]=node[c1].un[2];	    
      
    }
    //done setting up mymesh and mydata


    /* find the derivatives at the vertices by integrating over 8 subzonal blocks */
    
    /*block 1*/

    mydata.p[0]=0;
    mydata.p[1]=1;
    mydata.p[2]=4;
    mydata.p[3]=3;
    mydata.p[4]=9;
    mydata.p[5]=10;
    mydata.p[6]=13;
    mydata.p[7]=12;
    
    xc[0]=yc[0]=zc[0]=0.0;
    for (n=0;n<8;n++) {
      xc[0]+=0.125*mymesh.xc[mydata.p[n]];
      yc[0]+=0.125*mymesh.yc[mydata.p[n]];
      zc[0]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 2*/

    mydata.p[0]=1;
    mydata.p[1]=2;
    mydata.p[2]=5;
    mydata.p[3]=4;
    mydata.p[4]=10;
    mydata.p[5]=11;
    mydata.p[6]=14;
    mydata.p[7]=13;

    xc[1]=yc[1]=zc[1]=0.0;
    for (n=0;n<8;n++) {
      xc[1]+=0.125*mymesh.xc[mydata.p[n]];
      yc[1]+=0.125*mymesh.yc[mydata.p[n]];
      zc[1]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 3*/

    mydata.p[0]=10;
    mydata.p[1]=11;
    mydata.p[2]=14;
    mydata.p[3]=13;
    mydata.p[4]=19;
    mydata.p[5]=20;
    mydata.p[6]=23;
    mydata.p[7]=22;

    xc[2]=yc[2]=zc[2]=0.0;
    for (n=0;n<8;n++) {
      xc[2]+=0.125*mymesh.xc[mydata.p[n]];
      yc[2]+=0.125*mymesh.yc[mydata.p[n]];
      zc[2]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 4*/

    mydata.p[0]=9;
    mydata.p[1]=10;
    mydata.p[2]=13;
    mydata.p[3]=12;
    mydata.p[4]=18;
    mydata.p[5]=19;
    mydata.p[6]=22;
    mydata.p[7]=21;

    xc[3]=yc[3]=zc[3]=0.0;
    for (n=0;n<8;n++) {
      xc[3]+=0.125*mymesh.xc[mydata.p[n]];
      yc[3]+=0.125*mymesh.yc[mydata.p[n]];
      zc[3]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 5*/

    mydata.p[0]=3;
    mydata.p[1]=4;
    mydata.p[2]=7;
    mydata.p[3]=6;
    mydata.p[4]=12;
    mydata.p[5]=13;
    mydata.p[6]=16;
    mydata.p[7]=15;

    xc[4]=yc[4]=zc[4]=0.0;
    for (n=0;n<8;n++) {
      xc[4]+=0.125*mymesh.xc[mydata.p[n]];
      yc[4]+=0.125*mymesh.yc[mydata.p[n]];
      zc[4]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 6*/

    mydata.p[0]=4;
    mydata.p[1]=5;
    mydata.p[2]=8;
    mydata.p[3]=7;
    mydata.p[4]=13;
    mydata.p[5]=14;
    mydata.p[6]=17;
    mydata.p[7]=16;

    xc[5]=yc[5]=zc[5]=0.0;
    for (n=0;n<8;n++) {
      xc[5]+=0.125*mymesh.xc[mydata.p[n]];
      yc[5]+=0.125*mymesh.yc[mydata.p[n]];
      zc[5]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 7*/

    mydata.p[0]=13;
    mydata.p[1]=14;
    mydata.p[2]=17;
    mydata.p[3]=16;
    mydata.p[4]=22;
    mydata.p[5]=23;
    mydata.p[6]=26;
    mydata.p[7]=25;

    xc[6]=yc[6]=zc[6]=0.0;
    for (n=0;n<8;n++) {
      xc[6]+=0.125*mymesh.xc[mydata.p[n]];
      yc[6]+=0.125*mymesh.yc[mydata.p[n]];
      zc[6]+=0.125*mymesh.zc[mydata.p[n]];
    }

    block2(mymesh,mydata);

    /*block 8*/

    mydata.p[0]=12;
    mydata.p[1]=13;
    mydata.p[2]=16;
    mydata.p[3]=15;
    mydata.p[4]=21;
    mydata.p[5]=22;
    mydata.p[6]=25;
    mydata.p[7]=24;

    xc[7]=yc[7]=zc[7]=0.0;
    for (n=0;n<8;n++) {
      xc[7]+=0.125*mymesh.xc[mydata.p[n]];
      yc[7]+=0.125*mymesh.yc[mydata.p[n]];
      zc[7]+=0.125*mymesh.zc[mydata.p[n]];
    }
    
    block2(mymesh,mydata);

    /* find the derivates at the cell center */
  
    for (n=0;n<3;n++) {
      mydata.sx[n]/=mydata.vol;
      mydata.sy[n]/=mydata.vol;
      mydata.sz[n]/=mydata.vol;
    }
  

    /* loop over all 26 neighboring cells to find min and max values */

    fmax[0]=fmin[0]=mydata.u[13];
    fmax[1]=fmin[1]=mydata.v[13];
    fmax[2]=fmin[2]=mydata.w[13];

    for (n=0;n<27;n++) {
      
      if (mydata.u[n] > fmax[0])
	fmax[0]=mydata.u[n];
      else if (mydata.u[n] < fmin[0])
	fmin[0]=mydata.u[n];

      if (mydata.v[n] > fmax[1]) 
	fmax[1]=mydata.v[n];
      else if (mydata.v[n] < fmin[1])
	fmin[1]=mydata.v[n];
      
      if (mydata.w[n] > fmax[2]) 
	fmax[2]=mydata.w[n];
      else if (mydata.w[n] < fmin[2])
	fmin[2]=mydata.w[n];
    }

    fmax[0]-=mydata.u[13];
    fmax[1]-=mydata.v[13];
    fmax[2]-=mydata.w[13];

    fmin[0]-=mydata.u[13];
    fmin[1]-=mydata.v[13];
    fmin[2]-=mydata.w[13];

    /* apply Barth Jesperson limiter to the cell centered derivates by looping over 8 verticies l=1 to 8*/
    for (n=0;n<3;n++) {

      phi[n]=mylimiter(fmax[n],fmin[n],mydata.sx[n]*(xc[0]-mymesh.xc[13])
		       +mydata.sy[n]*(yc[0]-mymesh.yc[13])
		       +mydata.sz[n]*(zc[0]-mymesh.zc[13]));

      for (l=1;l<=7;l++) {

	factor=mylimiter(fmax[n],fmin[n],mydata.sx[n]*(xc[l]-mymesh.xc[13])
			 +mydata.sy[n]*(yc[l]-mymesh.yc[13])
			 +mydata.sz[n]*(zc[l]-mymesh.zc[13]));

	if (factor<phi[n]) phi[n]=factor;
      }
    }

    for (n=0;n<3;n++) {

      mydata.sx[n]*=phi[n];
      mydata.sy[n]*=phi[n];
      mydata.sz[n]*=phi[n];
      
    }

    /*u velocity gradient at vertex v*/
    node[v].sx_u[0]=mydata.sx[0];
    node[v].sx_u[1]=mydata.sy[0];
    node[v].sx_u[2]=mydata.sz[0];

    /*v velocity gradient at vertex v*/
    node[v].sx_v[0]=mydata.sx[1];
    node[v].sx_v[1]=mydata.sy[1];
    node[v].sx_v[2]=mydata.sz[1];

    /*w velocity gradient at vertex v*/
    node[v].sx_w[0]=mydata.sx[2];
    node[v].sx_w[1]=mydata.sy[2];
    node[v].sx_w[2]=mydata.sz[2];

   }// end if vflag equals 1

  } /* end vertex loop */

  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=72;
   
  int sendcounts[numtasks];
  int sdispls[numtasks];

  int total_counts=0;
   
  int num_reqs=0;
   
   kk=0;
   
   for (int j=0;j<numtasks;j++) { //loop over all ranks
     
     if (j!=rank ) {
       if (mcon.send_num[kk]>0) {
	 sendcounts[j]=ndat*mcon.send_num[kk];
	 sdispls[j]=total_counts;
	 total_counts+=ndat*mcon.send_num[kk];
	 num_reqs++;
       }
       else {
	 sendcounts[j]=0;
	 sdispls[j]=0;
       }
       kk++;       
     }
     else {
       sendcounts[j]=0;
       sdispls[j]=0;
     }
   }

   //create array holding data to send to other ranks
   send= new double[total_counts];

   //now fill the send array with the appropriate data
   //kk iterates over all other ranks in sequential order
   //k iterates over all cells on the current rank whose data is send to other ranks
   kk=k=0;
   for (int j=0;j<numtasks;j++) { //loop over all ranks
     
     int offset=sdispls[j];
     
     if (j!=rank) {

       if (mcon.send_num[kk]>0) {
	 
	 //copy local cell data into send array
	 for (int i=0;i<mcon.send_num[kk];i++,k++) {

	   auto c=cs[mcon.send_local[k]];
	   
	   int index=0;

	   //loop over all corners in the cell
	   for (auto cr : mesh.corners_of_cell(c)) {
	     
	     auto v = mesh.corners[cr].point;

	     //MUFx
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].sx_u[ii];
	       index++;

	     }
	     //MVFx
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].sx_v[ii];
	       index++;

	     }
	     //MWFx
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].sx_w[ii];
	       index++;

	     }	     

	   }//end corner loop	 


	 } //end cell loop
       } //end if mcon.send_num[kk]>0

       kk++;
       
     }//end if rank!=j
    
   } //done loop over ranks   
    
   //construct array to hold data being received from other ranks
   int recvcounts[numtasks];
   int rdispls[numtasks];
   
   kk2=total_counts=0;
   for (int j=0;j<numtasks;j++) { //loop over all ranks less than my rank

     if (rank!=j) {
       
       if (mcon.recv_num[kk2]>0) {
	 
	 recvcounts[j]=ndat*mcon.recv_num[kk2];
	 rdispls[j]=total_counts;
	 total_counts+=ndat*mcon.recv_num[kk2];
	 num_reqs++;
       }
       else {
	 recvcounts[j]=0;
	 rdispls[j]=0;
       }
       kk2++;       
     }
     else {
       recvcounts[j]=0;
       rdispls[j]=0;
     }
   }

   //create array holding data received from other ranks
   recv= new double[total_counts];


   MPI_Request *reqs;

   reqs= new  MPI_Request[num_reqs];

   //send requests
   kk=k=0;
   for (int j=0;j<numtasks;j++) { //loop over all ranks
     
     int offset=sdispls[j];
     int count=sendcounts[j];
     
     if (j!=rank) {

       if (count>0) {

	 MPI_Isend(send+offset, count, MPI_DOUBLE, j, 9, MPI_COMM_WORLD, &reqs[k]);

	 k++;
	 
       } //end if mcon.send_num[kk]>0
       
     }//end if rank!=j
    
   } //done loop over ranks

   //receive requests
   kk=0;
   for (int j=0;j<numtasks;j++) { //loop over all ranks
     
     int offset=rdispls[j];
     int count=recvcounts[j];
     
     if (j!=rank) {

       if (count>0) {

	 MPI_Irecv(recv+offset, count, MPI_DOUBLE, j, 9, MPI_COMM_WORLD, &reqs[k]);

	 k++;
	 
       } //end if mcon.send_num[kk]>0
       
     }//end if rank!=j
    
   } //done loop over ranks

   
   MPI_Waitall(k, reqs, MPI_STATUSES_IGNORE);

   delete reqs;

   //kk2 iterates over all other ranks in sequential order
   //k2 iterates over all cells on other ranks whose data is send to the current rank   
   kk2=k2=0;
   for (int j=0;j<numtasks;j++) { //loop over all ranks less than my rank
     
     int offset=rdispls[j];

     if (rank!=j) {
       
        if (mcon.recv_num[kk2]>0) {
	  
	 //copy data received from other rank into local ghost cells
	 for (int i=0;i<mcon.recv_num[kk2];i++,k2++) {

	   auto c=cs[mcon.recv_local[k2]];

	   int index=0;
	   
	   //find the correct order of vertices for cell c
	   myvert[0]=cell[c].iord1;
	   myvert[1]=cell[c].iord2;
	   myvert[2]=cell[c].iord3;
	   myvert[3]=cell[c].iord4;
	   myvert[4]=cell[c].iord5;
	   myvert[5]=cell[c].iord6;
	   myvert[6]=cell[c].iord7;
	   myvert[7]=cell[c].iord8;	

	   jj=0;
	   for ( auto v : mesh.vertices_of_cell(c) ) {

	     for (ii=0;ii<8;ii++)
	       if (node[v].gid==myvert[ii]) iorder[jj]=ii;

	     jj++;
	   }

	   ii=0;
	   for ( auto v : mesh.vertices_of_cell(c) ) {
	     cn[iorder[ii]]=v;
	     ii++;
	   }	   

	   //loop over corners 
	   for (int cr=0; cr<8; cr++) {
	     
	     auto v = verts[cn[cr]];

	     //MUFx
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].sx_u[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     //MVFx
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].sx_v[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     //MWFx
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].sx_w[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     
	   }//end corner loop	 


	 }//end cell loop

       } //end if recv_num[kk]>0
	 
       
       kk2++;

     } //end if rank!=j

   } //done loop over ranks


   delete[] send;
   delete[] recv;
   delete[] vflag;
  
}
