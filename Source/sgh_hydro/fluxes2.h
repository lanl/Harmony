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
void fluxes2( 
  mesh_t  &mesh,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {

  int i,j,k,l,n,kk,ii,jj;

  int imin,imax,jmin,jmax,kmin,kmax;

  int n1,n0,i1,i2,i3,i4,p;

  double vol,vol1;
  double dx,dy,dz;
  double x,y,z,mflux;

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;


  int iorder[8],cn[8],cn2[8],myvert[8],counter;

  int g1,g2,g3,g4;
  int h1,h2,h3,h4;
  int e1,e2;

  double mfac;

  const double TOL=1.e-25;

  int indx[27],mycells[8];

  //cells
  auto cs = mesh.cells();
  auto verts = mesh.vertices();
  auto NP = verts.size();
  auto cs_owned = mesh.cells_owned();
  
  int rank,numtasks,kk2,k2;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  
  double *recv,*send;

  MPI_Status stat;    

  mfac=1.0;

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

  //loop over all vertices 
  for ( auto v : verts ) {

    if (vflag[v]==1) { //skip ghost vertices

    //begin construction of the indx array containing vertex neighbors
    i=v;
    for (j=0;j<27;j++)
      indx[j]=-1; //default setting of -1 indicates boundary vertex

    indx[13]=v; //central vertex on the 27 point stencil

    for (j=0;j<8;j++)
      mycells[j]=-1; //default setting of -1 indicates boundary cell

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

        mycells[6]=c;

	indx[13]=cn[0];
	indx[14]=cn[1];
	indx[17]=cn[2];
	indx[16]=cn[3];
	indx[22]=cn[4];
	indx[23]=cn[5];
	indx[26]=cn[6];
	indx[25]=cn[7];


      }
      else if (l==1) { 

        mycells[7]=c;

	indx[12]=cn[0];
	indx[13]=cn[1];
	indx[16]=cn[2];
	indx[15]=cn[3];
	indx[21]=cn[4];
	indx[22]=cn[5];
	indx[25]=cn[6];
	indx[24]=cn[7];


      }

      else if (l==2) { 

        mycells[4]=c;

	indx[9]=cn[0];
	indx[10]=cn[1];
	indx[13]=cn[2];
	indx[12]=cn[3];
	indx[18]=cn[4];
	indx[19]=cn[5];
	indx[22]=cn[6];
	indx[21]=cn[7];


      }
      else if (l==3) { 

        mycells[5]=c;

	indx[10]=cn[0];
	indx[11]=cn[1];
	indx[14]=cn[2];
	indx[13]=cn[3];
	indx[19]=cn[4];
	indx[20]=cn[5];
	indx[23]=cn[6];
	indx[22]=cn[7];


      }

      else if (l==4) { 

        mycells[2]=c;

	indx[4]=cn[0];
	indx[5]=cn[1];
	indx[8]=cn[2];
	indx[7]=cn[3];
	indx[13]=cn[4];
	indx[14]=cn[5];
	indx[17]=cn[6];
	indx[16]=cn[7];


      }
      else if (l==5) { 

        mycells[3]=c;

	indx[3]=cn[0];
	indx[4]=cn[1];
	indx[7]=cn[2];
	indx[6]=cn[3];
	indx[12]=cn[4];
	indx[13]=cn[5];
	indx[16]=cn[6];
	indx[15]=cn[7];


      }
      else if (l==6) { 

        mycells[0]=c;

	indx[0]=cn[0];
	indx[1]=cn[1];
	indx[4]=cn[2];
	indx[3]=cn[3];
	indx[9]=cn[4];
	indx[10]=cn[5];
	indx[13]=cn[6];
	indx[12]=cn[7];


      }
      else if (l==7) { 

        mycells[1]=c;

	indx[1]=cn[0];
	indx[2]=cn[1];
	indx[5]=cn[2];
	indx[4]=cn[3];
	indx[10]=cn[4];
	indx[11]=cn[5];
	indx[14]=cn[6];
	indx[13]=cn[7];

      }

    } //finished looping over all cells containing vertex v


    if (indx[14]!=-1) { //compute u flux 

      /* vertex indicies */
      n0=indx[13];
      n1=indx[14];    
 
      /* indicies for cell centers which form verticies of the momentum control volume*/
      i1=mycells[2];    
      i2=mycells[6];   
      i3=mycells[5];   
      i4=mycells[1];


      /* xcp,ycp,zcp are the cell center positions before the Lagrange step */
      /* xc,yc,zc are the cell center positions after the Lagrange step */
      if (i1!=-1 && i2!=-1 && i3!=-1 && i4!=-1 ) { //interior vertex, not on a boundary
      
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if ( i1==-1 && i2==-1 && i3!=-1 && i4!=-1 ) { //vertex is on +x boundary

	//indices for first boundary face
        h1=indx[4];
	h2=indx[13];
	h3=indx[14];
	h4=indx[5];

	//indices for second boundary face
        g1=indx[13];
	g2=indx[22];
	g3=indx[23];
	g4=indx[14];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 2
	x2=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y2=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z2=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y6=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z6=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1!=-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on -x boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[22];
	h3=indx[23];
	h4=indx[14];
	
	//indices for second boundary face
        g1=indx[4];
	g2=indx[13];
	g3=indx[14];
	g4=indx[5];	
	
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
        // point 3
	x3=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y3=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z3=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y7=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z7=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }

       else if ( i1!=-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on +y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[17];
	h4=indx[14];

	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[14];
	g4=indx[11];

        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];	

        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1==-1 && i2!=-1 && i3!=-1 && i4==-1 ) { //vertex is on -y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[17];
	h4=indx[14];

	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[14];
	g4=indx[11];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];	

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }     
      else if (i1==-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[4];
	h2=indx[13];
	h3=indx[14];
	h4=indx[5];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[14];
	
	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[14];
	g4=indx[11];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y2=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z2=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];	
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y6=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z6=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	// point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];		

      }     
      else if (i1!=-1 && i2==-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[17];
	h4=indx[14];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[14];
	
	//indices for second boundary face
        g1=indx[4];
	g2=indx[13];
	g3=indx[14];
	g4=indx[5];
	
	// point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 3
	x3=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y3=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z3=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y7=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z7=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }     
      else if (i1==-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[17];
	h4=indx[14];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[14];
	
	//indices for second boundary face
        g1=indx[13];
	g2=indx[22];
	g3=indx[23];
	g4=indx[14];
	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

	// point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];

        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);	

        // point 4
	x4=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y4=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z4=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

        // point 8
	x8=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y8=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z8=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	

      }           
      else if (i1==-1 && i2==-1 && i3!=-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[22];
	h3=indx[23];
	h4=indx[14];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[14];
	
	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[14];
	g4=indx[11];
	

        // point 1
	x1=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y1=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z1=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	

	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);


        // point 5
	x5=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y5=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z5=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }

      /* flux  volume */

      vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
               x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
               x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
               x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
               x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
               x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
               x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
	    x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;

      /* centroid of the flux volume */

      x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
      y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
      z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);


      /* upwind test */
      if (vol<0.) { 
	p=n1;
	vol=-vol;
      }
      else {
	p=n0;
      }

      /* calculate the fluxes */
      if (vol>TOL) {

	/* x,y,z are the vertex position coordinates after the Lagrange step */

	dx=x-mesh.points[verts[p]].px[0];
	dy=y-mesh.points[verts[p]].px[1];
	dz=z-mesh.points[verts[p]].px[2];

	//compute mass flux for vertex v by taking 1/8th of the mass flux from each cell connected to v
	mflux=0.0;
	for (n=0;n<8;n++) {
	  if (mycells[n]!=-1) {
	    auto c1 = cs[mycells[n]];
	    mflux+=0.125*cell[c1].MFx[0];
	  }
	}


	   
	node[v].MUFx[0]=mflux*(node[verts[p]].un[0]+mfac*(dx*node[verts[p]].sx_u[0]+dy*node[verts[p]].sx_u[1]+dz*node[verts[p]].sx_u[2])); /* u momentum flux */
	node[v].MVFx[0]=mflux*(node[verts[p]].un[1]+mfac*(dx*node[verts[p]].sx_v[0]+dy*node[verts[p]].sx_v[1]+dz*node[verts[p]].sx_v[2])); /* v momentum flux */
	node[v].MWFx[0]=mflux*(node[verts[p]].un[2]+mfac*(dx*node[verts[p]].sx_w[0]+dy*node[verts[p]].sx_w[1]+dz*node[verts[p]].sx_w[2])); /* w momentum flux */

	//        if (rank==0) printf("v= %d flux= %e %e %e %e\n",(int)v,node[v].MUFx[0],node[v].MVFx[0],node[v].MWFx[0],mflux);

      }
      else {

	node[v].MUFx[0]=0.0; /* u momentum flux */
	node[v].MVFx[0]=0.0; /* v momentum flux */
	node[v].MWFx[0]=0.0; /* w momentum flux */

      }

    } //finished computing u fluxes

    
    if (indx[16]!=-1) { //compute v flux 

      /* vertex indicies */
      n0=indx[13];
      n1=indx[16];    
 
      /* indicies for cell centers which form verticies of the momentum control volume*/
      i1=mycells[3];    
      i2=mycells[7];    
      i3=mycells[6];    
      i4=mycells[2];    


      /* xcp,ycp,zcp are the cell center positions before the Lagrange step */
      /* xc,yc,zc are the cell center positions after the Lagrange step */
      if (i1!=-1 && i2!=-1 && i3!=-1 && i4!=-1 ) { //interior vertex, not on a boundary
      
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if ( i1==-1 && i2==-1 && i3!=-1 && i4!=-1 ) { //vertex is on +x boundary

	//indices for first boundary face
        h1=indx[4];
	h2=indx[13];
	h3=indx[16];
	h4=indx[7];

	//indices for second boundary face
        g1=indx[13];
	g2=indx[22];
	g3=indx[25];
	g4=indx[16];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 2
	x2=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y2=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z2=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y6=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z6=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1!=-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on -x boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[22];
	h3=indx[25];
	h4=indx[16];
	
	//indices for second boundary face
        g1=indx[4];
	g2=indx[13];
	g3=indx[16];
	g4=indx[7];	
	
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
        // point 3
	x3=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y3=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z3=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y7=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z7=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }

       else if ( i1!=-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on +y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[12];
	h3=indx[15];
	h4=indx[16];

	//indices for second boundary face
        g1=indx[14];
	g2=indx[13];
	g3=indx[16];
	g4=indx[17];

        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];	

        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1==-1 && i2!=-1 && i3!=-1 && i4==-1 ) { //vertex is on -y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[12];
	h3=indx[15];
	h4=indx[16];

	//indices for second boundary face
        g1=indx[14];
	g2=indx[13];
	g3=indx[16];
	g4=indx[17];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];	

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }     
      else if (i1==-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[4];
	h2=indx[13];
	h3=indx[16];
	h4=indx[7];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[16];
	
	//indices for second boundary face
        g1=indx[14];
	g2=indx[13];
	g3=indx[16];
	g4=indx[17];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y2=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z2=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];	
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y6=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z6=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	// point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];		

      }     
      else if (i1!=-1 && i2==-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[12];
	h3=indx[15];
	h4=indx[16];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[16];
	
	//indices for second boundary face
        g1=indx[4];
	g2=indx[13];
	g3=indx[16];
	g4=indx[7];
	
	// point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 3
	x3=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y3=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z3=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y7=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z7=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }     
      else if (i1==-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[12];
	h3=indx[15];
	h4=indx[16];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[16];
	
	//indices for second boundary face
        g1=indx[13];
	g2=indx[22];
	g3=indx[25];
	g4=indx[16];
	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

	// point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];

        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);	

        // point 4
	x4=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y4=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z4=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

        // point 8
	x8=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y8=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z8=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	

      }           
      else if (i1==-1 && i2==-1 && i3!=-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[22];
	h3=indx[25];
	h4=indx[16];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[16];
	
	//indices for second boundary face
        g1=indx[14];
	g2=indx[13];
	g3=indx[16];
	g4=indx[17];
	

        // point 1
	x1=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y1=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z1=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	

	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);


        // point 5
	x5=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y5=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z5=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }

      /* flux  volume */

      vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
               x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
               x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
               x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
               x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
               x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
               x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
	    x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;

      /* centroid of the flux volume */

      x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
      y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
      z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);


      /* upwind test */
      if (vol<0.) { 
	p=n1;
	vol=-vol;
       }
      else {
	p=n0;
      }

      /* calculate the fluxes */
      if (vol>TOL) {

	/* x,y,z are the vertex position coordinates after the Lagrange step */

	dx=x-mesh.points[verts[p]].px[0];
	dy=y-mesh.points[verts[p]].px[1];
	dz=z-mesh.points[verts[p]].px[2];

	//compute mass flux for vertex v by taking 1/8th of the mass flux from each cell connected to v
	mflux=0.0;
	for (n=0;n<8;n++) {
	  if (mycells[n]!=-1) {
	    auto c1 = cs[mycells[n]];
	    mflux+=0.125*cell[c1].MFx[1];
	  }
	}

	   
	node[v].MUFx[1]=mflux*(node[verts[p]].un[0]+mfac*(dx*node[verts[p]].sx_u[0]+dy*node[verts[p]].sx_u[1]+dz*node[verts[p]].sx_u[2])); /* u momentum flux */
	node[v].MVFx[1]=mflux*(node[verts[p]].un[1]+mfac*(dx*node[verts[p]].sx_v[0]+dy*node[verts[p]].sx_v[1]+dz*node[verts[p]].sx_v[2])); /* v momentum flux */
	node[v].MWFx[1]=mflux*(node[verts[p]].un[2]+mfac*(dx*node[verts[p]].sx_w[0]+dy*node[verts[p]].sx_w[1]+dz*node[verts[p]].sx_w[2])); /* w momentum flux */

      }
      else {

	node[v].MUFx[1]=0.0; /* u momentum flux */
	node[v].MVFx[1]=0.0; /* v momentum flux */
	node[v].MWFx[1]=0.0; /* w momentum flux */

      }

    } //finished computing v fluxes

    if (indx[22]!=-1) { //compute w flux 

      /* vertex indicies */
      n0=indx[13];
      n1=indx[22];    
 
      /* indicies for cell centers which form verticies of the momentum control volume*/
      i1=mycells[5];      /*cell i+1,j,k+1*/
      i2=mycells[6];      /*cell i+1,j+1,k+1*/
      i3=mycells[7];      /*cell i,j+1,k+1*/ 
      i4=mycells[4];      /*cell i,j,k+1*/ 


      /* xcp,ycp,zcp are the cell center positions before the Lagrange step */
      /* xc,yc,zc are the cell center positions after the Lagrange step */
      if (i1!=-1 && i2!=-1 && i3!=-1 && i4!=-1 ) { //interior vertex, not on a boundary
      
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if ( i1==-1 && i2==-1 && i3!=-1 && i4!=-1 ) { //vertex is on +x boundary

	//indices for first boundary face
        h1=indx[10];
	h2=indx[13];
	h3=indx[22];
	h4=indx[19];

	//indices for second boundary face
        g1=indx[13];
	g2=indx[16];
	g3=indx[25];
	g4=indx[22];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 2
	x2=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y2=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z2=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];

        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y6=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z6=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];

        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1!=-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on -x boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[25];
	h4=indx[22];
	
	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[22];
	g4=indx[19];	
	
        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
        // point 3
	x3=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y3=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z3=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y7=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z7=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }

       else if ( i1!=-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on +y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[14];
	h3=indx[23];
	h4=indx[22];

	//indices for second boundary face
        g1=indx[12];
	g2=indx[13];
	g3=indx[22];
	g4=indx[21];

        // point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
        // point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];

        // point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];	

        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);
	
        // point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];

      }
      else if (i1==-1 && i2!=-1 && i3!=-1 && i4==-1 ) { //vertex is on -y boundary

	//indices for first boundary face
        h1=indx[13];
	h2=indx[14];
	h3=indx[23];
	h4=indx[22];

	//indices for second boundary face
        g1=indx[12];
	g2=indx[13];
	g3=indx[22];
	g4=indx[21];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];	

        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

      }     
      else if (i1==-1 && i2==-1 && i3==-1 && i4!=-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[10];
	h2=indx[13];
	h3=indx[22];
	h4=indx[19];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[22];
	
	//indices for second boundary face
        g1=indx[12];
	g2=indx[13];
	g3=indx[22];
	g4=indx[21];	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 2
	x2=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y2=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z2=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 4
	x4=cell[cs[i4]].xcp[0];
	y4=cell[cs[i4]].xcp[1];
	z4=cell[cs[i4]].xcp[2];	
	
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 6
	x6=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y6=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z6=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	// point 8
	x8=cell[cs[i4]].xc[0];
	y8=cell[cs[i4]].xc[1];
	z8=cell[cs[i4]].xc[2];		

      }     
      else if (i1!=-1 && i2==-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[14];
	h3=indx[23];
	h4=indx[22];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[22];
	
	//indices for second boundary face
        g1=indx[10];
	g2=indx[13];
	g3=indx[22];
	g4=indx[19];
	
	// point 1
	x1=cell[cs[i1]].xcp[0];
	y1=cell[cs[i1]].xcp[1];
	z1=cell[cs[i1]].xcp[2];
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

        // point 3
	x3=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y3=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z3=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);

	// point 5
	x5=cell[cs[i1]].xc[0];
	y5=cell[cs[i1]].xc[1];
	z5=cell[cs[i1]].xc[2];
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);

        // point 7
	x7=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y7=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z7=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }     
      else if (i1==-1 && i2!=-1 && i3==-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[14];
	h3=indx[23];
	h4=indx[22];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[22];
	
	//indices for second boundary face
        g1=indx[13];
	g2=indx[16];
	g3=indx[25];
	g4=indx[22];
	
	
        // point 1
	x1=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y1=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z1=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);

	// point 2
	x2=cell[cs[i2]].xcp[0];
	y2=cell[cs[i2]].xcp[1];
	z2=cell[cs[i2]].xcp[2];

        // point 3
	x3=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y3=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z3=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);	

        // point 4
	x4=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y4=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z4=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);
	
        // point 5
	x5=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y5=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z5=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 6
	x6=cell[cs[i2]].xc[0];
	y6=cell[cs[i2]].xc[1];
	z6=cell[cs[i2]].xc[2];
	
        // point 7
	x7=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y7=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z7=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

        // point 8
	x8=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y8=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z8=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	

      }           
      else if (i1==-1 && i2==-1 && i3!=-1 && i4==-1 ) { //vertex is on a boundary edge 

	//indices for first boundary face
        h1=indx[13];
	h2=indx[16];
	h3=indx[25];
	h4=indx[22];

        //indices for first boundary edge
	e1=indx[13];
	e2=indx[22];
	
	//indices for second boundary face
        g1=indx[12];
	g2=indx[13];
	g3=indx[22];
	g4=indx[21];
	

        // point 1
	x1=0.5*(node[verts[e1]].xp[0]+node[verts[e2]].xp[0]);
	y1=0.5*(node[verts[e1]].xp[1]+node[verts[e2]].xp[1]);
	z1=0.5*(node[verts[e1]].xp[2]+node[verts[e2]].xp[2]);	
	
        // point 2
	x2=0.25*(node[verts[h1]].xp[0]+node[verts[h2]].xp[0]+
		 node[verts[h3]].xp[0]+node[verts[h4]].xp[0]);
	y2=0.25*(node[verts[h1]].xp[1]+node[verts[h2]].xp[1]+
		 node[verts[h3]].xp[1]+node[verts[h4]].xp[1]);
	z2=0.25*(node[verts[h1]].xp[2]+node[verts[h2]].xp[2]+
		 node[verts[h3]].xp[2]+node[verts[h4]].xp[2]);
	
	// point 3
	x3=cell[cs[i3]].xcp[0];
	y3=cell[cs[i3]].xcp[1];
	z3=cell[cs[i3]].xcp[2];	

	
        // point 4
	x4=0.25*(node[verts[g1]].xp[0]+node[verts[g2]].xp[0]+
		 node[verts[g3]].xp[0]+node[verts[g4]].xp[0]);
	y4=0.25*(node[verts[g1]].xp[1]+node[verts[g2]].xp[1]+
		 node[verts[g3]].xp[1]+node[verts[g4]].xp[1]);
	z4=0.25*(node[verts[g1]].xp[2]+node[verts[g2]].xp[2]+
		 node[verts[g3]].xp[2]+node[verts[g4]].xp[2]);


        // point 5
	x5=0.5*(mesh.points[verts[e1]].px[0]+mesh.points[verts[e2]].px[0]);
	y5=0.5*(mesh.points[verts[e1]].px[1]+mesh.points[verts[e2]].px[1]);
	z5=0.5*(mesh.points[verts[e1]].px[2]+mesh.points[verts[e2]].px[2]);	
	
        // point 6
	x6=0.25*(mesh.points[verts[h1]].px[0]+mesh.points[verts[h2]].px[0]+
		 mesh.points[verts[h3]].px[0]+mesh.points[verts[h4]].px[0]);
	y6=0.25*(mesh.points[verts[h1]].px[1]+mesh.points[verts[h2]].px[1]+
		 mesh.points[verts[h3]].px[1]+mesh.points[verts[h4]].px[1]);
	z6=0.25*(mesh.points[verts[h1]].px[2]+mesh.points[verts[h2]].px[2]+
		 mesh.points[verts[h3]].px[2]+mesh.points[verts[h4]].px[2]);
	
	// point 7
	x7=cell[cs[i3]].xc[0];
	y7=cell[cs[i3]].xc[1];
	z7=cell[cs[i3]].xc[2];
	
        // point 8
	x8=0.25*(mesh.points[verts[g1]].px[0]+mesh.points[verts[g2]].px[0]+
		 mesh.points[verts[g3]].px[0]+mesh.points[verts[g4]].px[0]);
	y8=0.25*(mesh.points[verts[g1]].px[1]+mesh.points[verts[g2]].px[1]+
		 mesh.points[verts[g3]].px[1]+mesh.points[verts[g4]].px[1]);
	z8=0.25*(mesh.points[verts[g1]].px[2]+mesh.points[verts[g2]].px[2]+
		 mesh.points[verts[g3]].px[2]+mesh.points[verts[g4]].px[2]);

	

      }
      
      /* flux  volume */
      
	vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
               x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
               x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
               x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
               x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
               x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
               x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
	    x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;

      /* centroid of the flux volume */

      x=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
      y=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
      z=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);


      /* upwind test */
      if (vol<0.) { 
	p=n1;
	vol=-vol;
      }
      else {
	p=n0;
      }

      /* calculate the fluxes */
      if (vol>TOL) {

	/* x,y,z are the vertex position coordinates after the Lagrange step */

	dx=x-mesh.points[verts[p]].px[0];
	dy=y-mesh.points[verts[p]].px[1];
	dz=z-mesh.points[verts[p]].px[2];

	//compute mass flux for vertex v by taking 1/8th of the mass flux from each cell connected to v
	mflux=0.0;
	for (n=0;n<8;n++) {
	  if (mycells[n]!=-1) {
	    auto c1 = cs[mycells[n]];
	    mflux+=0.125*cell[c1].MFx[2];
	  }
	}

	   
	node[v].MUFx[2]=mflux*(node[verts[p]].un[0]+mfac*(dx*node[verts[p]].sx_u[0]+dy*node[verts[p]].sx_u[1]+dz*node[verts[p]].sx_u[2])); /* u momentum flux */
	node[v].MVFx[2]=mflux*(node[verts[p]].un[1]+mfac*(dx*node[verts[p]].sx_v[0]+dy*node[verts[p]].sx_v[1]+dz*node[verts[p]].sx_v[2])); /* v momentum flux */
	node[v].MWFx[2]=mflux*(node[verts[p]].un[2]+mfac*(dx*node[verts[p]].sx_w[0]+dy*node[verts[p]].sx_w[1]+dz*node[verts[p]].sx_w[2])); /* w momentum flux */

      }
      else {

	node[v].MUFx[2]=0.0; /* u momentum flux */
	node[v].MVFx[2]=0.0; /* v momentum flux */
	node[v].MWFx[2]=0.0; /* w momentum flux */

      }

    } //finished computing w fluxes

    } //end vflag equals one

  } //done vertex loop


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
	      
	       send[ndat*i+index+offset]=node[v].MUFx[ii];
	       index++;

	     }
	     //MVFx
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].MVFx[ii];
	       index++;

	     }
	     //MWFx
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].MWFx[ii];
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
	       if (vflag[v]==0) node[v].MUFx[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     //MVFx
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].MVFx[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     //MWFx
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].MWFx[ii]=recv[ndat*i+index+offset];
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
