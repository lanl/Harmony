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
class unstruct {

public :

  int nregs;
  int kkzl;
  int kkzll;
  int kksl;
  int kksll;
  int kkmatlistll;
  int* kkz_zs_start;
  int* kkzs_s;
  int* kstyp;
  int* kkss2;
  int* kksz;
  int* kksz2;  
  int* kksp1;
  int* kksp2;
  int* kksf;
  int* kkse;
  int* ktagsz;
  int* nummatsz;

  int kkpll;
  int kkfll;
  int* kkfz;
  int* kztyp;
  int kkell;
  int* ketyp;

  double* zzv;
  double* zx;
  double* fx;
  double* px;
  double* ex;

  double* dpx;
  double* dex;
  double* dfx;
  double* dzx;
  double* fdv;
  double* sdvz;

  int* kkz_to_mpyg;

};

double tripleproduct(double* x,double* y,double* z)
{
  
  double tp;

  tp = 0.0;
  tp = tp + x[0]*(y[1]*z[2]-y[2]*z[1]);
  tp = tp + x[1]*(y[2]*z[0]-y[0]*z[2]);
  tp = tp + x[2]*(y[0]*z[1]-y[1]*z[0]);

  return tp;
}

void fluxes1( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {

  int i,j,k,n,indx;

  int imin,imax,jmin,jmax,kmin,kmax;

  int n1,n0,i1,i2,i3,i4,p,cz;

  double vol,vol1;
  double dx,dy,dz;
  double x,y,z;

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  int kk;
  int hpyg,mpyg,mpygc;
  int imatcnt,maxpriority;

  int p1,p2,f,ii;
  int nstart,n2,np;

  int iorder[8],cn[8],cn2[8],myvert[8],counter;

  int h,s,s2;

  int ff,e,zz,zz2;
  double vv1,vv2,vv3,sdeltav,sdeltav1,sdeltav2;
  
  double v1a1b[3],v2a2b[3],vfafb[3],vzazb[3],veaeb[3];
  double vfa2a[3],vfa1a[3],v1afb[3],v1a2a[3],v2afb[3],veaza[3];
  double veafb[3],vfaea[3],vfaza[3],vzafb[3],vzaeb[3],vea2a[3],v2aeb[3];
  double v1aea[3],v1aeb[3];

  int sdonor,zdonor;
  int ipriority,imatcount;

  double mfac;


  const double TOL=1.e-10;
  const double TOLV=1.e-6;

  //cells
  auto cs = mesh.cells();
  auto verts = mesh.vertices();
  auto faces = mesh.list_faces();
  auto edges = mesh.list_edges();

  auto myface = faces[0];

  auto NC=cs.size();
  auto NP=verts.size();
  auto NF=faces.size();
  auto NE=edges.size();

  unstruct m;

  int mflag;

  double normal;

  int rank,numtasks,kk2,k2;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  
  double *recv,*send;

  MPI_Status stat;  

  mfac=1.0;
  
  m.kkzll=m.kkzl=NC;
  m.kkpll=NP;
  m.kksll=m.kksl=24*NC;
  m.kkfll=NF;
  m.kkell=NE;

  //side data structures
  m.kstyp=new int[m.kksll];
  m.kkss2=new int[m.kksll];
  m.kksz=new int[m.kksll];
  m.kksz2=new int[m.kksll];
  m.kksp1=new int[m.kksll];
  m.kksp2=new int[m.kksll];
  m.kksf=new int[m.kksll];
  m.kkse=new int[m.kksll];
  m.ktagsz=new int[m.kksll];
  
  m.sdvz=new double[m.kksll];

  //cell data structures
  m.kztyp=new int[m.kkzll];
  m.kkz_zs_start=new int[m.kkzll];
  m.nummatsz=new int[m.kkzll];
  m.kkz_to_mpyg=new int[m.kkzll];
  m.zzv=new double[m.kkzll];
  m.zx=new double[3*m.kkzll];
  m.dzx=new double[3*m.kkzll];

  //node data structures  
  m.px=new double[3*m.kkpll];
  m.dpx=new double[3*m.kkpll];

  //face data structures
  m.kkfz=new int[m.kkfll];
  m.fx=new double[3*m.kkfll];
  m.dfx=new double[3*m.kkfll];

  //edge data structures
  m.ex=new double[3*m.kkell];
  m.dex=new double[3*m.kkell];


  
  //loop over all nodes
  i=0;
  for ( auto v : mesh.vertices() ) {

    m.px[3*i]=mesh.points[v].px[0];
    m.px[1+3*i]=mesh.points[v].px[1];
    m.px[2+3*i]=mesh.points[v].px[2];

    m.dpx[3*i]=node[v].xp[0]-mesh.points[v].px[0];
    m.dpx[1+3*i]=node[v].xp[1]-mesh.points[v].px[1];
    m.dpx[2+3*i]=node[v].xp[2]-mesh.points[v].px[2];

    i++;
  }

  //loop over all faces
  for ( auto f : mesh.list_faces() ) {

    h=f;
    m.fx[3*h]=m.fx[1+3*h]=m.fx[2+3*h]=0.0;
    m.dfx[3*h]=m.dfx[1+3*h]=m.dfx[2+3*h]=0.0;
    for ( auto v : mesh.points_of_face(f) ) {

      for (j=0;j<3;j++) {
	m.fx[j+3*h]+=0.25*mesh.points[v].px[j];
	m.dfx[j+3*h]+=0.25*(node[v].xp[j]-mesh.points[v].px[j]);
      }
	  
	 
    }     
  }//end face loop

  //loop over all edges
  for ( auto f : mesh.list_edges() ) {

    h=f;
    m.ex[3*h]=m.ex[1+3*h]=m.ex[2+3*h]=0.0;
    m.dex[3*h]=m.dex[1+3*h]=m.dex[2+3*h]=0.0;
    for ( auto v : mesh.points_of_edge(f) ) {

      for (j=0;j<3;j++) {
	m.ex[j+3*h]+=0.5*mesh.points[v].px[j];
	m.dex[j+3*h]+=0.5*(node[v].xp[j]-mesh.points[v].px[j]);
      }
	  
	 
    }     
  }//end edge loop  

  /* loop over all cells to setup side data structures */
  s=0;
  for ( auto c : mesh.cells() ) {

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

        //construct array of cell vertices, cn, with vertices in proper order
        //find centroids for the cell, zx,dzx
	i=0;
	h=c;
	m.kztyp[h]=1;
	m.zx[3*h]=m.zx[1+3*h]=m.zx[2+3*h]=0.0;
	m.dzx[3*h]=m.dzx[1+3*h]=m.dzx[2+3*h]=0.0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
	  cn[iorder[i]]=v;

	  for (j=0;j<3;j++) {
	    m.zx[j+3*h]+=0.125*mesh.points[v].px[j];
	    m.dzx[j+3*h]+=0.125*(node[v].xp[j]-mesh.points[v].px[j]);
	  }
	  
	  i++;
	}
	
	i=c;
	m.kkz_zs_start[i]=s; //store the index of the first side for the cell 
	m.zzv[i]=cell[c].vol; //store cell volume
	
	//find +x face
	i1=myvert[0]=cn[5];
	i2=myvert[1]=cn[6];
	i3=myvert[2]=cn[2];
	i4=myvert[3]=cn[1];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;
	
	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[5];
	m.kksp2[s]=myvert[1]=cn[6];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[6];
	m.kksp2[s]=myvert[1]=cn[2];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[2];
	m.kksp2[s]=myvert[1]=cn[1];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[1];
	m.kksp2[s]=myvert[1]=cn[5];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//find +y face
	i1=myvert[0]=cn[6];
	i2=myvert[1]=cn[7];
	i3=myvert[2]=cn[3];
	i4=myvert[3]=cn[2];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;
	
	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[6];
	m.kksp2[s]=myvert[1]=cn[7];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[7];
	m.kksp2[s]=myvert[1]=cn[3];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[3];
	m.kksp2[s]=myvert[1]=cn[2];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[2];
	m.kksp2[s]=myvert[1]=cn[6];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;

	//find -x face
	i1=myvert[0]=cn[7];
	i2=myvert[1]=cn[4];
	i3=myvert[2]=cn[0];
	i4=myvert[3]=cn[3];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;

	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[7];
	m.kksp2[s]=myvert[1]=cn[4];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[4];
	m.kksp2[s]=myvert[1]=cn[0];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[0];
	m.kksp2[s]=myvert[1]=cn[3];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[3];
	m.kksp2[s]=myvert[1]=cn[7];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;	


	//find -y face
	i1=myvert[0]=cn[4];
	i2=myvert[1]=cn[5];
	i3=myvert[2]=cn[1];
	i4=myvert[3]=cn[0];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;
	
	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[4];
	m.kksp2[s]=myvert[1]=cn[5];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[5];
	m.kksp2[s]=myvert[1]=cn[1];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[1];
	m.kksp2[s]=myvert[1]=cn[0];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[0];
	m.kksp2[s]=myvert[1]=cn[4];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;	
	    
	//find +z face
	i1=myvert[0]=cn[7];
	i2=myvert[1]=cn[6];
	i3=myvert[2]=cn[5];
	i4=myvert[3]=cn[4];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;
	
	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[7];
	m.kksp2[s]=myvert[1]=cn[6];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[6];
	m.kksp2[s]=myvert[1]=cn[5];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[5];
	m.kksp2[s]=myvert[1]=cn[4];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[4];
	m.kksp2[s]=myvert[1]=cn[7];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//find -z face
	i1=myvert[0]=cn[0];
	i2=myvert[1]=cn[1];
	i3=myvert[2]=cn[2];
	i4=myvert[3]=cn[3];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}
	
	m.kkfz[myface]=c;
	
	//loop over all edges in the face
	//edge 1
	m.kksp1[s]=myvert[0]=cn[0];
	m.kksp2[s]=myvert[1]=cn[1];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 2
	m.kksp1[s]=myvert[0]=cn[1];
	m.kksp2[s]=myvert[1]=cn[2];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 3
	m.kksp1[s]=myvert[0]=cn[2];
	m.kksp2[s]=myvert[1]=cn[3];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;
	
	//edge 4
	m.kksp1[s]=myvert[0]=cn[3];
	m.kksp2[s]=myvert[1]=cn[0];
        //find edge associated with these points, store index in ee
	for (auto f : mesh.edges_of_face(myface)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_edge(f))
	    for (j=0;j<2;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==2) m.kkse[s]=f;
	}
	
	if (indx!=-1) {
	  m.kksz2[s]=indx;
	  m.kstyp[s]=1;
	}
	else {
	  m.kksz2[s]=-1;
	  m.kstyp[s]=0;
	}
	m.ktagsz[s]=m.kstyp[s];

	m.kksz[s]=c;
	m.kksf[s]=myface;

	s=s+1;

  } //end cell loop and setup of side data structures

  //setup kkss2 data structure
  for (s=0;s<m.kksll;s++) {
    m.kkss2[s]=-1;

    h=m.kksz2[s];
    if (h>-1) {
      //loop over all sides in zone h
      for (i =  m.kkz_zs_start[h]; i<(m.kkz_zs_start[h]+24);i++) 
	if ((m.kkse[s]==m.kkse[i]) && (m.kksf[s]==m.kksf[i])) m.kkss2[s]=i;
    }


  }//end setup of kkss2


  //calculate side advection flux volumes
  //       sdvz(s) > 0 => a volume sdvz(s) is advected
  //                      from zone kksz(s) into zone kksz2(s).
  //       sdvz(s) < 0 => a volume |sdvz(s)| is advected
  //                      from zone kksz2(s) into zone kksz(s).
  for (s=0; s<m.kksl; s++) {
    m.sdvz[s]=0.0;
  }
 
  for (s=0; s<m.kksl; s++) {
    if (m.kstyp[s]>0) {
      e = m.kkse[s];
      ff = m.kksf[s];
      zz = m.kksz[s];
 
      if (m.kkfz[ff]==zz) {
	p1 = m.kksp1[s];
	p2 = m.kksp2[s];
	s2 = m.kkss2[s];

	//Construct tetrahedron edge vectors
	for (j=0; j<3; j++) {
	  v1a1b[j] = m.dpx[j+3*p1];	
	  v2a2b[j] = m.dpx[j+3*p2];
	  veaeb[j] = m.dex[j+3*e];
	  vfafb[j] = m.dfx[j+3*ff];
	  vfa2a[j] = m.px[j+3*p2]-m.fx[j+3*ff];
	  vfa1a[j] = m.px[j+3*p1]-m.fx[j+3*ff];
	  v1afb[j] = m.fx[j+3*ff]+m.dfx[j+3*ff]-m.px[j+3*p1];
	  v1a2a[j] = m.px[j+3*p2]-m.px[j+3*p1];
	  v2afb[j] = m.fx[j+3*ff] + m.dfx[j+3*ff] - m.px[j+3*p2];
	  vfaea[j] = m.ex[j+3*e] - m.fx[j+3*ff];
	  veafb[j] = m.fx[j+3*ff] + m.dfx[j+3*ff]  - m.ex[j+3*e];
	  v1aea[j] = m.ex[j+3*e]-m.px[j+3*p1];
	  v1aeb[j] = m.ex[j+3*e] + m.dex[j+3*e] - m.px[j+3*p1];
	  vea2a[j] = m.px[j+3*p2]-m.ex[j+3*e];
	  v2aeb[j] = m.ex[j+3*e] + m.dex[j+3*e] - m.px[j+3*p2];
	}

	vv1 = tripleproduct(vfafb,vfaea,vfa2a);
	vv2 = tripleproduct(veaeb,vea2a,veafb);
	vv3 = tripleproduct(v2a2b,v2afb,v2aeb);
	sdeltav2 = (vv1 + vv2 + vv3)/6.0;

	vv1 = tripleproduct(vfafb,vfa1a,vfaea);
	vv2 = tripleproduct(v1aeb,v1aea,v1afb);
	vv3 = tripleproduct(v1aeb,v1afb,v1a1b);
	sdeltav1 = (vv1 + vv2 + vv3)/6.0;
	sdeltav  = sdeltav1 + sdeltav2;

	if (sdeltav>0.0) {
	  m.sdvz[s]   = sdeltav;
	  m.sdvz[s2]  = 0.0;
	}
	else {
	  m.sdvz[s]   = 0.0;
	  m.sdvz[s2]  = -sdeltav;
	}

      } //end if kkfz[f]==z
    } //end if m.m.kstyp[s]!=0 && m.m.ktagsz[s]>0
  } //end side loop
  
  
  /* loop over all cells to find the fluxes MFx, EFx, TFx */
  for ( auto c : mesh.cells_owned() ) {

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

        //construct array of cell vertices, cn, with vertices in proper order
	i=0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
	  cn[iorder[i]]=v;
	  i++;
	}


	//find +x face
	i1=myvert[0]=cn[2];
	i2=myvert[1]=cn[6];
	i3=myvert[2]=cn[5];
	i4=myvert[3]=cn[1];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}

        if (indx!=-1) { //construct x fluxes for interior face
 
          /* cell indicies on either side of the face */
	  n0 = c;
	  n1 = indx;
	  
	  /* xp,yp,zp are the vertex positions before the Lagrange step */
	  /* x,y,z are the vertex positions after the Lagrange step */
	  x1=node[i1].xp[0];
	  x2=node[i2].xp[0];
	  x3=node[i3].xp[0];
	  x4=node[i4].xp[0];
	  x5=mesh.points[i1].px[0];
	  x6=mesh.points[i2].px[0];
	  x7=mesh.points[i3].px[0];
	  x8=mesh.points[i4].px[0];

	  y1=node[i1].xp[1];
	  y2=node[i2].xp[1];
	  y3=node[i3].xp[1];
	  y4=node[i4].xp[1];
	  y5=mesh.points[i1].px[1];
	  y6=mesh.points[i2].px[1];
	  y7=mesh.points[i3].px[1];
	  y8=mesh.points[i4].px[1];

	  z1=node[i1].xp[2];
	  z2=node[i2].xp[2];
	  z3=node[i3].xp[2];
	  z4=node[i4].xp[2];
	  z5=mesh.points[i1].px[2];
	  z6=mesh.points[i2].px[2];
	  z7=mesh.points[i3].px[2];
	  z8=mesh.points[i4].px[2];
	  
	  /* centroid of the flux volume */
	  x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
	  y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
	  z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);

	  /* loop over all sides in the cell to find those sides corresponding to myface */
	  /* and add contributions from those sides to flux volume for myface */
	  vol=0.0;
	  for (s=m.kkz_zs_start[c]; s<(m.kkz_zs_start[c]+24);s++) {

	    s2=m.kkss2[s];

	    if (s2==-1) continue;

	    if (m.kksf[s]==myface && m.sdvz[s]>0.0) {
	      vol+=m.sdvz[s]; //region flows out of the cell
	    }
	    else if (m.kksf[s2]==myface && m.sdvz[s2]>0.0) {  
	      vol-=m.sdvz[s2]; //region flows into the cell
	    }

	  }//end loop over sides in the cell	  


	  /* upwind test */
	  if (vol<0.) {
	    p=n1;
            vol1=-vol;
	  }
	  else {
	    p=n0;
            vol1=vol;	    
	  }

	  cz = cs[p];
	  dx=x-cell[cz].xc[0];
	  dy=y-cell[cz].xc[1];
	  dz=z-cell[cz].xc[2];

	  cell[c].MFx[0]=vol*(cell[cz].dc+mfac*(dx*cell[cz].sx_d[0]+dy*cell[cz].sx_d[1]+dz*cell[cz].sx_d[2])); /* mass */
	  cell[c].EFx[0]=cell[c].MFx[0]*(cell[cz].ec+mfac*(dx*cell[cz].sx_ie[0]+dy*cell[cz].sx_ie[1]+dz*cell[cz].sx_ie[2])); /* internal energy */
	  cell[c].TFx[0]=cell[c].MFx[0]*(cell[cz].enc+mfac*(dx*cell[cz].sx_en[0]+dy*cell[cz].sx_en[1]+dz*cell[cz].sx_en[2])); /* total energy */	  
	  

	} // if indx not equal to -1
        // end construction of flux on +x face


	//find +y face
	i1=myvert[0]=cn[3];
	i2=myvert[1]=cn[7];
	i3=myvert[2]=cn[6];
	i4=myvert[3]=cn[2];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}

    
        if (indx!=-1) { //construct y flux for interior face
 
          /* cell indicies on either side of the face */
	  n0 = c;
	  n1 = indx;

	  /* xp,yp,zp are the vertex positions before the Lagrange step */
	  /* x,y,z are the vertex positions after the Lagrange step */	  
	  x1=node[i1].xp[0];
	  x2=node[i2].xp[0];
	  x3=node[i3].xp[0];
	  x4=node[i4].xp[0];
	  x5=mesh.points[i1].px[0];
	  x6=mesh.points[i2].px[0];
	  x7=mesh.points[i3].px[0];
	  x8=mesh.points[i4].px[0];

	  y1=node[i1].xp[1];
	  y2=node[i2].xp[1];
	  y3=node[i3].xp[1];
	  y4=node[i4].xp[1];
	  y5=mesh.points[i1].px[1];
	  y6=mesh.points[i2].px[1];
	  y7=mesh.points[i3].px[1];
	  y8=mesh.points[i4].px[1];

	  z1=node[i1].xp[2];
	  z2=node[i2].xp[2];
	  z3=node[i3].xp[2];
	  z4=node[i4].xp[2];
	  z5=mesh.points[i1].px[2];
	  z6=mesh.points[i2].px[2];
	  z7=mesh.points[i3].px[2];
	  z8=mesh.points[i4].px[2];

	  /* centroid of the flux volume */
	  x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
	  y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
	  z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);

	  /* loop over all sides in the cell to find those sides corresponding to myface */
	  /* and add contributions from those sides to flux volume for myface */
	  vol=0.0;
	  for (s=m.kkz_zs_start[c]; s<(m.kkz_zs_start[c]+24);s++) {

	    s2=m.kkss2[s];

	    if (s2==-1) continue;

	    if (m.kksf[s]==myface && m.sdvz[s]>0.0) {
	      vol+=m.sdvz[s]; //region flows out of the cell
	    }
	    else if (m.kksf[s2]==myface && m.sdvz[s2]>0.0) {  
	      vol-=m.sdvz[s2]; //region flows into the cell
	    }

	  }//end loop over sides in the cell	  


	  /* upwind test */
	  if (vol<0.) {
	    p=n1;
            vol1=-vol;
	  }
	  else {
	    p=n0;
            vol1=vol;	    
	  }

	  cz = cs[p];
	  dx=x-cell[cz].xc[0];
	  dy=y-cell[cz].xc[1];
	  dz=z-cell[cz].xc[2];

	  cell[c].MFx[1]=vol*(cell[cz].dc+mfac*(dx*cell[cz].sx_d[0]+dy*cell[cz].sx_d[1]+dz*cell[cz].sx_d[2])); /* mass */
	  cell[c].EFx[1]=cell[c].MFx[1]*(cell[cz].ec+mfac*(dx*cell[cz].sx_ie[0]+dy*cell[cz].sx_ie[1]+dz*cell[cz].sx_ie[2])); /* internal energy */
	  cell[c].TFx[1]=cell[c].MFx[1]*(cell[cz].enc+mfac*(dx*cell[cz].sx_en[0]+dy*cell[cz].sx_en[1]+dz*cell[cz].sx_en[2])); /* total energy */	  
	  	  

	} // if indx not equal to -1
        // end construction of flux on +y face


	//find +z face
	i1=myvert[0]=cn[5];
	i2=myvert[1]=cn[6];
	i3=myvert[2]=cn[7];
	i4=myvert[3]=cn[4];

        indx=-1;
        i=c;
	for (auto f : mesh.faces_of_cell(c)) {
	  counter=0;
	  for (auto v1 : mesh.points_of_face(f))
	    for (j=0;j<4;j++)
	      if (myvert[j]==v1) counter++;
	  if (counter==4) myface=f;
	}

	for (auto c1 : mesh.cells_of_face(myface)) {
          j=c1;
	  if (j!=i) indx=j;  //if indx is not equal to -1 then face is not on boundary
	}

    
        if (indx!=-1) { //construct y flux for interior face
 
          /* cell indicies on either side of the face */
	  n0 = c;
	  n1 = indx;

	  /* xp,yp,zp are the vertex positions before the Lagrange step */
	  /* x,y,z are the vertex positions after the Lagrange step */
	  x1=node[i1].xp[0];
	  x2=node[i2].xp[0];
	  x3=node[i3].xp[0];
	  x4=node[i4].xp[0];
	  x5=mesh.points[i1].px[0];
	  x6=mesh.points[i2].px[0];
	  x7=mesh.points[i3].px[0];
	  x8=mesh.points[i4].px[0];

	  y1=node[i1].xp[1];
	  y2=node[i2].xp[1];
	  y3=node[i3].xp[1];
	  y4=node[i4].xp[1];
	  y5=mesh.points[i1].px[1];
	  y6=mesh.points[i2].px[1];
	  y7=mesh.points[i3].px[1];
	  y8=mesh.points[i4].px[1];

	  z1=node[i1].xp[2];
	  z2=node[i2].xp[2];
	  z3=node[i3].xp[2];
	  z4=node[i4].xp[2];
	  z5=mesh.points[i1].px[2];
	  z6=mesh.points[i2].px[2];
	  z7=mesh.points[i3].px[2];
	  z8=mesh.points[i4].px[2];

	  /* centroid of the flux volume */
	  x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
	  y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
	  z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);

	  /* loop over all sides in the cell to find those sides corresponding to myface */
	  /* and add contributions from those sides to flux volume for myface */
	  vol=0.0;
	  for (s=m.kkz_zs_start[c]; s<(m.kkz_zs_start[c]+24);s++) {

	    s2=m.kkss2[s];

	    if (s2==-1) continue;

	    if (m.kksf[s]==myface && m.sdvz[s]>0.0) {
	      vol+=m.sdvz[s]; //region flows out of the cell
	    }
	    else if (m.kksf[s2]==myface && m.sdvz[s2]>0.0) {  
	      vol-=m.sdvz[s2]; //region flows into the cell
	    }

	  }//end loop over sides in the cell	  


	  /* upwind test */
	  if (vol<0.) {
	    p=n1;
            vol1=-vol;
	  }
	  else {
	    p=n0;
            vol1=vol;	    
	  }

	  cz = cs[p];
	  dx=x-cell[cz].xc[0];
	  dy=y-cell[cz].xc[1];
	  dz=z-cell[cz].xc[2];

	  cell[c].MFx[2]=vol*(cell[cz].dc+mfac*(dx*cell[cz].sx_d[0]+dy*cell[cz].sx_d[1]+dz*cell[cz].sx_d[2])); /* mass */
	  cell[c].EFx[2]=cell[c].MFx[2]*(cell[cz].ec+mfac*(dx*cell[cz].sx_ie[0]+dy*cell[cz].sx_ie[1]+dz*cell[cz].sx_ie[2])); /* internal energy */
	  cell[c].TFx[2]=cell[c].MFx[2]*(cell[cz].enc+mfac*(dx*cell[cz].sx_en[0]+dy*cell[cz].sx_en[1]+dz*cell[cz].sx_en[2])); /* total energy */	  
	  



	} // if indx not equal to -1
        // end construction of flux on +z face


  } //end cell loop


  //free memory
  //side data structures
  delete[] m.kstyp;
  delete[] m.kkss2;
  delete[] m.kksz;
  delete[] m.kksz2;
  delete[] m.kksp1;
  delete[] m.kksp2;
  delete[] m.kksf;
  delete[] m.kkse;
  delete[] m.ktagsz;
  delete[] m.sdvz;
  
  //cell data structures
  delete[] m.kztyp;
  delete[] m.kkz_zs_start;
  delete[] m.nummatsz;
  delete[] m.kkz_to_mpyg;
  delete[] m.zzv;
  delete[] m.zx;
  delete[] m.dzx;

  //node data structures  
  delete[] m.px;
  delete[] m.dpx;

  //face data structures
  delete[] m.kkfz;
  delete[] m.fx;
  delete[] m.dfx;

  //edge data structures
  delete[] m.ex;
  delete[] m.dex;


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
	     send[ndat*i+index+offset]=cell[c].MFx[n];
	     index++;
	     send[ndat*i+index+offset]=cell[c].EFx[n];
	     index++;
	     send[ndat*i+index+offset]=cell[c].TFx[n];
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
	     cell[c].MFx[n]=recv[ndat*i+index+offset];
	     index++;
	     cell[c].EFx[n]=recv[ndat*i+index+offset];
	     index++;
	     cell[c].TFx[n]=recv[ndat*i+index+offset];
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
