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
double time(
  mesh_t  &mesh,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  double DT0,
  double DT1  
) {


  int i,j,k,n;

  int k1,k2,j1,j2,i1,i2;

  double E,sum;

  double x[8],y[8],z[8];
  double u[8],v[8],w[8];

  double dx,dy,dz;

  double xa,ya,za;

  double ua,va,wa;

  double s,ds;
  double CFL=0.25;
  double DTM;

  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();
  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  double DT;

  int iorder[8],cn[8],myvert[8];
  double ave[3];

  DT=DT1;

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
	ua=va=wa=0.0;
	for ( auto v : mesh.vertices_of_cell(c) ) {
	  x[iorder[i]]=mesh.points[v].px[0];
	  y[iorder[i]]=mesh.points[v].px[1];
	  z[iorder[i]]=mesh.points[v].px[2];
	  i++;
	  ua+=node[v].un[0];
	  va+=node[v].un[1];
	  wa+=node[v].un[2];
	}

	/*calculate first diagonal */
	dx=(x[5]-x[3]);
	dy=(y[5]-y[3]);
	dz=(z[5]-z[3]);

	ds=sqrt(dx*dx+dy*dy+dz*dz);

	s=fabs(ua*dx+va*dy+wa*dz)/ds;

	DTM=CFL*ds/(s+cell[c].sound);
	if (DTM<DT) {
	  DT=DTM;
	}


	/*calculate second diagonal */
	dx=(x[6]-x[0]);
	dy=(y[6]-y[0]);
	dz=(z[6]-z[0]);

	ds=sqrt(dx*dx+dy*dy+dz*dz);

	s=fabs(ua*dx+va*dy+wa*dz)/ds;

	DTM=CFL*ds/(s+cell[c].sound);
	if (DTM<DT) {
	  DT=DTM;
	}


	/*calculate third diagonal */
	dx=(x[7]-x[1]);
	dy=(y[7]-y[1]);
	dz=(z[7]-z[1]);

	ds=sqrt(dx*dx+dy*dy+dz*dz);

	s=fabs(ua*dx+va*dy+wa*dz)/ds;

	DTM=CFL*ds/(s+cell[c].sound);
	if (DTM<DT) {
	  DT=DTM;
	}


	/*calculate fourth diagonal */
	dx=(x[4]-x[2]);
	dy=(y[4]-y[2]);
	dz=(z[4]-z[2]);

	ds=sqrt(dx*dx+dy*dy+dz*dz);

	s=fabs(ua*dx+va*dy+wa*dz)/ds;

	DTM=CFL*ds/(s+cell[c].sound);
	if (DTM<DT) {
	  DT=DTM;
	}



  }

  // printf("DT= %f\n",DT);
 DTM=1.1*DT0;

 if (DT>DTM) DT=DTM;
 
  double time_step=DT;

  return time_step;


}






