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
void ale_init2( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node	    
) {

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  int i,j,k,l,n;
  double ke,sum,sum1;

  int iorder[8],cn[8],myvert[8];

  const double TOL=1.e-6;

  //vertices
  auto verts = mesh.vertices();

  //cells
  auto cs = mesh.cells();


  //find cell volume after the Lagrangian step
  for ( auto c : cs ) {

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

        //find cell centroids
        cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;
        cell[c].xcp[0]=cell[c].xcp[1]=cell[c].xcp[2]=0.0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
          for (j=0;j<3;j++) {
	    cell[c].xc[j]+=0.125*mesh.points[v].px[j];
            cell[c].xcp[j]+=0.125*node[v].xp[j];
	  }
	}

	x1=mesh.points[cn[0]].px[0];
	x2=mesh.points[cn[1]].px[0];
	x3=mesh.points[cn[2]].px[0];
	x4=mesh.points[cn[3]].px[0];
	x5=mesh.points[cn[4]].px[0];
	x6=mesh.points[cn[5]].px[0];
	x7=mesh.points[cn[6]].px[0];
	x8=mesh.points[cn[7]].px[0];

	y1=mesh.points[cn[0]].px[1];
	y2=mesh.points[cn[1]].px[1];
	y3=mesh.points[cn[2]].px[1];
	y4=mesh.points[cn[3]].px[1];
	y5=mesh.points[cn[4]].px[1];
	y6=mesh.points[cn[5]].px[1];
	y7=mesh.points[cn[6]].px[1];
	y8=mesh.points[cn[7]].px[1];

	z1=mesh.points[cn[0]].px[2];
	z2=mesh.points[cn[1]].px[2];
	z3=mesh.points[cn[2]].px[2];
	z4=mesh.points[cn[3]].px[2];
	z5=mesh.points[cn[4]].px[2];
	z6=mesh.points[cn[5]].px[2];
	z7=mesh.points[cn[6]].px[2];
	z8=mesh.points[cn[7]].px[2];


	/* cell volume */

        cell[c].vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
                 x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
                 x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
                 x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
                 x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
                 x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
                 x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
		 x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;

	//update the density at cell center
	cell[c].dc=cell[c].mass/cell[c].vol;

        //find kinetic energy at cell center
        ke=0.0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
	  ke+=0.125*0.5*node[v].mp*(node[v].un[0]*node[v].un[0]+node[v].un[1]*node[v].un[1]+node[v].un[2]*node[v].un[2]);
	}

	//find total energy at the cell center
	cell[c].enc=ke/cell[c].mass+cell[c].ec;

	// zero all cell-centered fluxes and gradients 
	for (n=0;n<=2;n++) {
	  cell[c].MFx[n]=cell[c].EFx[n]=cell[c].TFx[n]=0.0;
	  cell[c].sx_d[n]=cell[c].sx_ie[n]=cell[c].sx_en[n]=0.0;
	}

  }




}



