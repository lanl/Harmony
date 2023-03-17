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
#include <map>

double mylimiter(double amax, double amin, double c)
{

  if (c>0)
    return MIN(1.,amax/c);
  else if (c<0)
    return MIN(1.,amin/c);
  else
    return 1.;

}

typedef struct _mydata { /* data type used for computing gradients in "slopes1" and "block1" */

  double vol,sx[3],sy[3],sz[3];
  int p[8], id;
  double d[27],ie[27],en[27];

} mydata_t;  


void block1(mymesh_t & mesh,mat_t & mat,mydata_t & mydata) /*subroutine used by "slopes1" to find cell-centered gradients*/
{
  int i,j,k;
  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  double b11,b12,b13,b14,b15,b16,b17,b18;
  double b21,b22,b23,b24,b25,b26,b27,b28;
  double b31,b32,b33,b34,b35,b36,b37,b38;

  /*centroids of the neighboring cells forming the block 1 control volume*/

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

  /* volume of block 1 */

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

  /* block 1 density derivative */
  mydata.sx[0]+=(mydata.d[mydata.p[0]]*b11+mydata.d[mydata.p[1]]*b12+mydata.d[mydata.p[2]]*b13+mydata.d[mydata.p[3]]*b14+mydata.d[mydata.p[4]]*b15+mydata.d[mydata.p[5]]*b16+mydata.d[mydata.p[6]]*b17+mydata.d[mydata.p[7]]*b18);
  mydata.sy[0]+=(mydata.d[mydata.p[0]]*b21+mydata.d[mydata.p[1]]*b22+mydata.d[mydata.p[2]]*b23+mydata.d[mydata.p[3]]*b24+mydata.d[mydata.p[4]]*b25+mydata.d[mydata.p[5]]*b26+mydata.d[mydata.p[6]]*b27+mydata.d[mydata.p[7]]*b28);
  mydata.sz[0]+=(mydata.d[mydata.p[0]]*b31+mydata.d[mydata.p[1]]*b32+mydata.d[mydata.p[2]]*b33+mydata.d[mydata.p[3]]*b34+mydata.d[mydata.p[4]]*b35+mydata.d[mydata.p[5]]*b36+mydata.d[mydata.p[6]]*b37+mydata.d[mydata.p[7]]*b38);

  /* block 1 internal energy derivative */
  mydata.sx[1]+=(mydata.ie[mydata.p[0]]*b11+mydata.ie[mydata.p[1]]*b12+mydata.ie[mydata.p[2]]*b13+mydata.ie[mydata.p[3]]*b14+mydata.ie[mydata.p[4]]*b15+mydata.ie[mydata.p[5]]*b16+mydata.ie[mydata.p[6]]*b17+mydata.ie[mydata.p[7]]*b18);
  mydata.sy[1]+=(mydata.ie[mydata.p[0]]*b21+mydata.ie[mydata.p[1]]*b22+mydata.ie[mydata.p[2]]*b23+mydata.ie[mydata.p[3]]*b24+mydata.ie[mydata.p[4]]*b25+mydata.ie[mydata.p[5]]*b26+mydata.ie[mydata.p[6]]*b27+mydata.ie[mydata.p[7]]*b28);
  mydata.sz[1]+=(mydata.ie[mydata.p[0]]*b31+mydata.ie[mydata.p[1]]*b32+mydata.ie[mydata.p[2]]*b33+mydata.ie[mydata.p[3]]*b34+mydata.ie[mydata.p[4]]*b35+mydata.ie[mydata.p[5]]*b36+mydata.ie[mydata.p[6]]*b37+mydata.ie[mydata.p[7]]*b38);

  /* block 1 total energy derivative */
  mydata.sx[2]+=(mydata.en[mydata.p[0]]*b11+mydata.en[mydata.p[1]]*b12+mydata.en[mydata.p[2]]*b13+mydata.en[mydata.p[3]]*b14+mydata.en[mydata.p[4]]*b15+mydata.en[mydata.p[5]]*b16+mydata.en[mydata.p[6]]*b17+mydata.en[mydata.p[7]]*b18);
  mydata.sy[2]+=(mydata.en[mydata.p[0]]*b21+mydata.en[mydata.p[1]]*b22+mydata.en[mydata.p[2]]*b23+mydata.en[mydata.p[3]]*b24+mydata.en[mydata.p[4]]*b25+mydata.en[mydata.p[5]]*b26+mydata.en[mydata.p[6]]*b27+mydata.en[mydata.p[7]]*b28);
  mydata.sz[2]+=(mydata.en[mydata.p[0]]*b31+mydata.en[mydata.p[1]]*b32+mydata.en[mydata.p[2]]*b33+mydata.en[mydata.p[3]]*b34+mydata.en[mydata.p[4]]*b35+mydata.en[mydata.p[5]]*b36+mydata.en[mydata.p[6]]*b37+mydata.en[mydata.p[7]]*b38);

}


void slopes1(
  mesh_t  &mesh,
  mat_t &mat,
  global_data_t &gd,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {	     

  /* finds mydata.cell centered derivates of density, internal energy and stress deviators */

  double factor,phi[4],fmin[4],fmax[4];

  int n,l,i,j,k;
  int flag;

  mydata_t mydata; // data structure containing gradients 
  mymesh_t mymesh; // local mesh containing the 27 cell stencil 

  int iorder[8],cn[8],cn2[8],myvert[8],counter;

  const double TOL=1.e-6;

  int rank,numtasks,kk,kk2,k2;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
 
  double *recv,*send;

  MPI_Status stat;

  //cells
  auto cs = mesh.cells();
  auto verts = mesh.vertices();
  auto faces = mesh.list_faces();

  auto c1 = cs[0];
  auto c2 = cs[0];
  auto c3 = cs[0];
  auto c4 = cs[0];

  auto myface = faces[0];

  int indx2[27];

  std::map<int,int> global_to_local;
  
  for ( auto c : cs ) 
    global_to_local[cell[c].id]=c;
    

  /* loop over all cells */
  for ( auto c : mesh.cells_owned()) {

        myvert[0]=cell[c].iord1;
        myvert[1]=cell[c].iord2;
        myvert[2]=cell[c].iord3;
        myvert[3]=cell[c].iord4;
        myvert[4]=cell[c].iord5;
        myvert[5]=cell[c].iord6;
        myvert[6]=cell[c].iord7;
        myvert[7]=cell[c].iord8;

	j=0;
	for ( auto v : mesh.vertices_of_cell(c) ) {

	  for (i=0;i<8;i++) 
            if (node[v].gid==myvert[i]) iorder[j]=i;

	  j++;
	}

    
	i=0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
	  cn[iorder[i]]=v;
	  i++;
	}

        //begin construction of the indx array containing cells neighbors
        i=c;
        for (j=0;j<27;j++)
          mymesh.indx[j]=-1; //default setting of -1 indicates boundary cell

        mymesh.indx[13]=c; //central cell on the 27 point stencil
	  
	int iz = cell[c].id % (gd.LIMX*gd.LIMY);

	int ix = iz % gd.LIMX;

	int iy = (iz - ix)/gd.LIMX;

	iz = (cell[c].id -iy*gd.LIMX -ix)/(gd.LIMX*gd.LIMY);

	ix++;
	iy++;  
	iz++;

        int il=0;
        for (int kk=iz-1; kk<=iz+1; kk++)	
	  for (int jj=iy-1; jj<=iy+1; jj++)	
	    for (int ii=ix-1; ii<=ix+1; ii++) {
	      
	      n=INDEX_FOR_3D(ii,1,gd.LIMX,jj,1,gd.LIMY,kk,1);
	      
	      if (ii<1 || ii>gd.LIMX || jj<1 || jj>gd.LIMY || kk<1 || kk>gd.LIMZ)
		mymesh.indx[il]=-1;
	      else
		mymesh.indx[il]=global_to_local[n];

	      //if (indx2[il]!=mymesh.indx[il]) printf("rank= %d cell= %d indx= %d %d\n",rank,c,mymesh.indx[il],indx2[il]);
	      il++;
	    }      
	//end construction of the indx array containing cells neighbors
	
        //fill in the centroid array
        //fill all cells on the stencil whose index is not -1
        for (j=0;j<27;j++) 
          if (mymesh.indx[j]!=-1) {
	    c2 = cs[mymesh.indx[j]];
	    mymesh.xc[j]=cell[c2].xc[0];
	    mymesh.yc[j]=cell[c2].xc[1];
	    mymesh.zc[j]=cell[c2].xc[2];
	  }

	//fill the ghost and boundary cells on the stencil
	complete_stencil(mymesh);
 
	//setup mydata
	mydata.vol=0.0;

	for (n=0;n<3;n++) {
	  mydata.sx[n]=0.0;
	  mydata.sy[n]=0.0;
	  mydata.sz[n]=0.0;
	}

	//copy region data onto local mesh
	for (n=0;n<27;n++) {

	  c1=cs[mymesh.cpy[n]];
	  mydata.d[n]=cell[c1].dc;
	  mydata.ie[n]=cell[c1].ec;
	  mydata.en[n]=cell[c1].enc;
	  	    
	}

	/* find the derivatives at the cell center by integrating over 8 subzonal blocks */

	/*block 1*/

	mydata.p[0]=0;
	mydata.p[1]=1;
	mydata.p[2]=4;
	mydata.p[3]=3;
	mydata.p[4]=9;
	mydata.p[5]=10;
	mydata.p[6]=13;
	mydata.p[7]=12;

	block1(mymesh,mat,mydata);

	/*block 2*/

	mydata.p[0]=1;
	mydata.p[1]=2;
	mydata.p[2]=5;
	mydata.p[3]=4;
	mydata.p[4]=10;
	mydata.p[5]=11;
	mydata.p[6]=14;
	mydata.p[7]=13;

	block1(mymesh,mat,mydata);

	/*block 3*/

	mydata.p[0]=10;
	mydata.p[1]=11;
	mydata.p[2]=14;
	mydata.p[3]=13;
	mydata.p[4]=19;
	mydata.p[5]=20;
	mydata.p[6]=23;
	mydata.p[7]=22;

	block1(mymesh,mat,mydata);

	/*block 4*/

	mydata.p[0]=9;
	mydata.p[1]=10;
	mydata.p[2]=13;
	mydata.p[3]=12;
	mydata.p[4]=18;
	mydata.p[5]=19;
	mydata.p[6]=22;
	mydata.p[7]=21;

	block1(mymesh,mat,mydata);

	/*block 5*/

	mydata.p[0]=3;
	mydata.p[1]=4;
	mydata.p[2]=7;
	mydata.p[3]=6;
	mydata.p[4]=12;
	mydata.p[5]=13;
	mydata.p[6]=16;
	mydata.p[7]=15;

	block1(mymesh,mat,mydata);

	/*block 6*/

	mydata.p[0]=4;
	mydata.p[1]=5;
	mydata.p[2]=8;
	mydata.p[3]=7;
	mydata.p[4]=13;
	mydata.p[5]=14;
	mydata.p[6]=17;
	mydata.p[7]=16;

	block1(mymesh,mat,mydata);

	/*block 7*/

	mydata.p[0]=13;
	mydata.p[1]=14;
	mydata.p[2]=17;
	mydata.p[3]=16;
	mydata.p[4]=22;
	mydata.p[5]=23;
	mydata.p[6]=26;
	mydata.p[7]=25;

	block1(mymesh,mat,mydata);

	/*block 8*/

	mydata.p[0]=12;
	mydata.p[1]=13;
	mydata.p[2]=16;
	mydata.p[3]=15;
	mydata.p[4]=21;
	mydata.p[5]=22;
	mydata.p[6]=25;
	mydata.p[7]=24;
	
	block1(mymesh,mat,mydata);

	/* find the derivates at the cell center */
  
	for (n=0;n<3;n++) {
	  mydata.sx[n]/=mydata.vol;
	  mydata.sy[n]/=mydata.vol;
	  mydata.sz[n]/=mydata.vol;
	}
  

	/* loop over all 26 neighboring cells to find min and max values */

	fmax[0]=fmin[0]=mydata.d[13];
	fmax[1]=fmin[1]=mydata.ie[13];
	fmax[2]=fmin[2]=mydata.en[13];
	
	for (n=0;n<27;n++) {

	  if (mydata.d[n] > fmax[0])
	    fmax[0]=mydata.d[n];
	  else if (mydata.d[n] < fmin[0])
	    fmin[0]=mydata.d[n];

	  if (mydata.ie[n] > fmax[1]) 
	    fmax[1]=mydata.ie[n];
	  else if (mydata.ie[n] < fmin[1])
	    fmin[1]=mydata.ie[n];

	  if (mydata.en[n] > fmax[2]) 
	    fmax[2]=mydata.en[n];
	  else if (mydata.en[n] < fmin[2])
	    fmin[2]=mydata.en[n];

	}

	fmax[0]-=mydata.d[13];
	fmax[1]-=mydata.ie[13];
	fmax[2]-=mydata.en[13];
	
	fmin[0]-=mydata.d[13];
	fmin[1]-=mydata.ie[13];
	fmin[2]-=mydata.en[13];
	
	/* apply Barth Jesperson limiter to the cell centered derivates by looping over 8 verticies l=1 to 8*/

	for (n=0;n<3;n++) {
	      
	  phi[n]=mylimiter(fmax[n],fmin[n],mydata.sx[n]*(mesh.points[cn[0]].px[0]-mymesh.xc[13])
			   +mydata.sy[n]*(mesh.points[cn[0]].px[1]-mymesh.yc[13])
			   +mydata.sz[n]*(mesh.points[cn[0]].px[2]-mymesh.zc[13]));

	  for (l=1;l<=7;l++) {

	    factor=mylimiter(fmax[n],fmin[n],mydata.sx[n]*(mesh.points[cn[l]].px[0]-mymesh.xc[13])
			     +mydata.sy[n]*(mesh.points[cn[l]].px[1]-mymesh.yc[13])
			     +mydata.sz[n]*(mesh.points[cn[l]].px[2]-mymesh.zc[13]));
	    
	    if (factor<phi[n]) phi[n]=factor;
	  }
	}

	for (n=0;n<3;n++) {

	  mydata.sx[n]*=phi[n];
	  mydata.sy[n]*=phi[n];
	  mydata.sz[n]*=phi[n];
	  
	}

	/*density gradient at cell center */
	cell[c].sx_d[0]=mydata.sx[0];
	cell[c].sx_d[1]=mydata.sy[0];
	cell[c].sx_d[2]=mydata.sz[0];

	/*internal energy gradient at cell center */
	cell[c].sx_ie[0]=mydata.sx[1];
	cell[c].sx_ie[1]=mydata.sy[1];
	cell[c].sx_ie[2]=mydata.sz[1];

	/*total energy gradient at cell center */
	cell[c].sx_en[0]=mydata.sx[2];
	cell[c].sx_en[1]=mydata.sy[2];
	cell[c].sx_en[2]=mydata.sz[2];

		
  } //end cell loop

  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=9;
   
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

	   for (int n=0;n<3;n++) {
	     send[ndat*i+index+offset]=cell[c].sx_d[n];
	     index++;
	     send[ndat*i+index+offset]=cell[c].sx_ie[n];
	     index++;
	     send[ndat*i+index+offset]=cell[c].sx_en[n];
	     index++;
	   }

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
	   
	   for (int n=0;n<3;n++) {
	     cell[c].sx_d[n]=recv[ndat*i+index+offset];
	     index++;
	     cell[c].sx_ie[n]=recv[ndat*i+index+offset];
	     index++;
	     cell[c].sx_en[n]=recv[ndat*i+index+offset];
	     index++;
	   }

	 }//end cell loop

       } //end if recv_num[kk]>0
	 
       
       kk2++;

     } //end if rank!=j

   } //done loop over ranks


   delete[] send;
   delete[] recv;


 

}
