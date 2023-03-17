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

void metrics( 
  mesh_t  &mesh,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node
) {

  int i,j,k,n;

  int k1,k2,j1,j2,i1,i2;

  double E,sum;

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;

  double u1,u2,u3,u4,u5,u6,u7,u8;
  double v1,v2,v3,v4,v5,v6,v7,v8;
  double w1,w2,w3,w4,w5,w6,w7,w8;

  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();
  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  double temp;
  
  int numtasks,rank;
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  int iorder[8],cn[8],myvert[8];
  double ave[3];

  double factor,tau10[3],tau20[3],tau30[3];

  /* calculate the cell volume, density and forces at the vertices */
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



    x1=mesh.points[verts[cn[0]]].px[0];
    x2=mesh.points[verts[cn[1]]].px[0];
    x3=mesh.points[verts[cn[2]]].px[0];
    x4=mesh.points[verts[cn[3]]].px[0];
    x5=mesh.points[verts[cn[4]]].px[0];
    x6=mesh.points[verts[cn[5]]].px[0];
    x7=mesh.points[verts[cn[6]]].px[0];
    x8=mesh.points[verts[cn[7]]].px[0];

    y1=mesh.points[verts[cn[0]]].px[1];
    y2=mesh.points[verts[cn[1]]].px[1];
    y3=mesh.points[verts[cn[2]]].px[1];
    y4=mesh.points[verts[cn[3]]].px[1];
    y5=mesh.points[verts[cn[4]]].px[1];
    y6=mesh.points[verts[cn[5]]].px[1];
    y7=mesh.points[verts[cn[6]]].px[1];
    y8=mesh.points[verts[cn[7]]].px[1];

    z1=mesh.points[verts[cn[0]]].px[2];
    z2=mesh.points[verts[cn[1]]].px[2];
    z3=mesh.points[verts[cn[2]]].px[2];
    z4=mesh.points[verts[cn[3]]].px[2];
    z5=mesh.points[verts[cn[4]]].px[2];
    z6=mesh.points[verts[cn[5]]].px[2];
    z7=mesh.points[verts[cn[6]]].px[2];
    z8=mesh.points[verts[cn[7]]].px[2];


    /* cell volume */
    cell[c].vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
                   x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) + 
                   x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
                   x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
                   x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
                   x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
                   x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
		x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;

    /* cell centroid */
    cell[c].xc[0]=0.125*(x1+x2+x3+x4+x5+x6+x7+x8);
    cell[c].xc[1]=0.125*(y1+y2+y3+y4+y5+y6+y7+y8);
    cell[c].xc[2]=0.125*(z1+z2+z3+z4+z5+z6+z7+z8);



    /* cell B matrix */
    cell[c].b1[0]=(-y2*z3 -y2*z4 +y2*z5 +y2*z6 +y3*z2 -y3*z4 
+y4*z2 +y4*z3 -y4*z5 -y4*z8 -y5*z2 +y5*z4 
-y5*z6 +y5*z8 -y6*z2 +y6*z5 +y8*z4 -y8*z5 )/12.;

    cell[c].b2[0]=(+y1*z3 +y1*z4 -y1*z5 -y1*z6 -y3*z1 -y3*z4 
+y3*z6 +y3*z7 -y4*z1 +y4*z3 +y5*z1 -y5*z6 
+y6*z1 -y6*z3 +y6*z5 -y6*z7 -y7*z3 +y7*z6 )/12.;

    cell[c].b3[0]=(-y1*z2 +y1*z4 +y2*z1 +y2*z4 -y2*z6 -y2*z7 
-y4*z1 -y4*z2 +y4*z7 +y4*z8 +y6*z2 -y6*z7 
+y7*z2 -y7*z4 +y7*z6 -y7*z8 -y8*z4 +y8*z7 )/12.;

    cell[c].b4[0]=(-y1*z2 -y1*z3 +y1*z5 +y1*z8 +y2*z1 -y2*z3 
+y3*z1 +y3*z2 -y3*z7 -y3*z8 -y5*z1 +y5*z8 
+y7*z3 -y7*z8 -y8*z1 +y8*z3 -y8*z5 +y8*z7 )/12.;

    cell[c].b5[0]=(+y1*z2 -y1*z4 +y1*z6 -y1*z8 -y2*z1 +y2*z6 
+y4*z1 -y4*z8 -y6*z1 -y6*z2 +y6*z7 +y6*z8 
-y7*z6 +y7*z8 +y8*z1 +y8*z4 -y8*z6 -y8*z7 )/12.;

    cell[c].b6[0]=(+y1*z2 -y1*z5 -y2*z1 +y2*z3 -y2*z5 +y2*z7 
-y3*z2 +y3*z7 +y5*z1 +y5*z2 -y5*z7 -y5*z8 
-y7*z2 -y7*z3 +y7*z5 +y7*z8 +y8*z5 -y8*z7)/12.;

    cell[c].b7[0]=(+y2*z3 -y2*z6 -y3*z2 +y3*z4 -y3*z6 +y3*z8 
-y4*z3 +y4*z8 +y5*z6 -y5*z8 +y6*z2 +y6*z3 
-y6*z5 -y6*z8 -y8*z3 -y8*z4 +y8*z5 +y8*z6 )/12.;

    cell[c].b8[0]=(-y1*z4 +y1*z5 +y3*z4 -y3*z7 +y4*z1 -y4*z3 
+y4*z5 -y4*z7 -y5*z1 -y5*z4 +y5*z6 +y5*z7 
-y6*z5 +y6*z7 +y7*z3 +y7*z4 -y7*z5 -y7*z6 )/12.;


    cell[c].b1[1]=(-z2*x3 -z2*x4 +z2*x5 +z2*x6 +z3*x2 -z3*x4 
+z4*x2 +z4*x3 -z4*x5 -z4*x8 -z5*x2 +z5*x4 
-z5*x6 +z5*x8 -z6*x2 +z6*x5 +z8*x4 -z8*x5 )/12.;

    cell[c].b2[1]=(+z1*x3 +z1*x4 -z1*x5 -z1*x6 -z3*x1 -z3*x4 
+z3*x6 +z3*x7 -z4*x1 +z4*x3 +z5*x1 -z5*x6 
+z6*x1 -z6*x3 +z6*x5 -z6*x7 -z7*x3 +z7*x6 )/12.;

    cell[c].b3[1]=(-z1*x2 +z1*x4 +z2*x1 +z2*x4 -z2*x6 -z2*x7 
-z4*x1 -z4*x2 +z4*x7 +z4*x8 +z6*x2 -z6*x7 
+z7*x2 -z7*x4 +z7*x6 -z7*x8 -z8*x4 +z8*x7 )/12.;

    cell[c].b4[1]=(-z1*x2 -z1*x3 +z1*x5 +z1*x8 +z2*x1 -z2*x3 
+z3*x1 +z3*x2 -z3*x7 -z3*x8 -z5*x1 +z5*x8 
+z7*x3 -z7*x8 -z8*x1 +z8*x3 -z8*x5 +z8*x7 )/12.;

    cell[c].b5[1]=(+z1*x2 -z1*x4 +z1*x6 -z1*x8 -z2*x1 +z2*x6 
+z4*x1 -z4*x8 -z6*x1 -z6*x2 +z6*x7 +z6*x8 
-z7*x6 +z7*x8 +z8*x1 +z8*x4 -z8*x6 -z8*x7 )/12.;

    cell[c].b6[1]=(+z1*x2 -z1*x5 -z2*x1 +z2*x3 -z2*x5 +z2*x7 
-z3*x2 +z3*x7 +z5*x1 +z5*x2 -z5*x7 -z5*x8 
-z7*x2 -z7*x3 +z7*x5 +z7*x8 +z8*x5 -z8*x7 )/12.;

    cell[c].b7[1]=(+z2*x3 -z2*x6 -z3*x2 +z3*x4 -z3*x6 +z3*x8 
-z4*x3 +z4*x8 +z5*x6 -z5*x8 +z6*x2 +z6*x3 
-z6*x5 -z6*x8 -z8*x3 -z8*x4 +z8*x5 +z8*x6 )/12.;

    cell[c].b8[1]=(-z1*x4 +z1*x5 +z3*x4 -z3*x7 +z4*x1 -z4*x3 
+z4*x5 -z4*x7 -z5*x1 -z5*x4 +z5*x6 +z5*x7 
-z6*x5 +z6*x7 +z7*x3 +z7*x4 -z7*x5 -z7*x6 )/12.;


    cell[c].b1[2]=(-x2*y3 -x2*y4 +x2*y5 +x2*y6 +x3*y2 -x3*y4 
+x4*y2 +x4*y3 -x4*y5 -x4*y8 -x5*y2 +x5*y4 
-x5*y6 +x5*y8 -x6*y2 +x6*y5 +x8*y4 -x8*y5 )/12.;

    cell[c].b2[2]=(+x1*y3 +x1*y4 -x1*y5 -x1*y6 -x3*y1 -x3*y4 
+x3*y6 +x3*y7 -x4*y1 +x4*y3 +x5*y1 -x5*y6 
+x6*y1 -x6*y3 +x6*y5 -x6*y7 -x7*y3 +x7*y6 )/12.;

    cell[c].b3[2]=(-x1*y2 +x1*y4 +x2*y1 +x2*y4 -x2*y6 -x2*y7 
-x4*y1 -x4*y2 +x4*y7 +x4*y8 +x6*y2 -x6*y7 
+x7*y2 -x7*y4 +x7*y6 -x7*y8 -x8*y4 +x8*y7 )/12.;

    cell[c].b4[2]=(-x1*y2 -x1*y3 +x1*y5 +x1*y8 +x2*y1 -x2*y3 
+x3*y1 +x3*y2 -x3*y7 -x3*y8 -x5*y1 +x5*y8 
+x7*y3 -x7*y8 -x8*y1 +x8*y3 -x8*y5 +x8*y7 )/12.;

    cell[c].b5[2]=(+x1*y2 -x1*y4 +x1*y6 -x1*y8 -x2*y1 +x2*y6 
+x4*y1 -x4*y8 -x6*y1 -x6*y2 +x6*y7 +x6*y8 
-x7*y6 +x7*y8 +x8*y1 +x8*y4 -x8*y6 -x8*y7 )/12.;

    cell[c].b6[2]=(+x1*y2 -x1*y5 -x2*y1 +x2*y3 -x2*y5 +x2*y7 
-x3*y2 +x3*y7 +x5*y1 +x5*y2 -x5*y7 -x5*y8 
-x7*y2 -x7*y3 +x7*y5 +x7*y8 +x8*y5 -x8*y7 )/12.;

    cell[c].b7[2]=(+x2*y3 -x2*y6 -x3*y2 +x3*y4 -x3*y6 +x3*y8 
-x4*y3 +x4*y8 +x5*y6 -x5*y8 +x6*y2 +x6*y3 
-x6*y5 -x6*y8 -x8*y3 -x8*y4 +x8*y5 +x8*y6 )/12.;

    cell[c].b8[2]=(-x1*y4 +x1*y5 +x3*y4 -x3*y7 +x4*y1 -x4*y3 
+x4*y5 -x4*y7 -x5*y1 -x5*y4 +x5*y6 +x5*y7 
-x6*y5 +x6*y7 +x7*y3 +x7*y4 -x7*y5 -x7*y6 )/12.;



    /* compute the velocity gradient */
    for (j=0;j<3;j++) {

      cell[c].dudx1[j]=(node[verts[cn[0]]].un[0]*cell[c].b1[j]+node[verts[cn[1]]].un[0]*cell[c].b2[j]+node[verts[cn[2]]].un[0]*cell[c].b3[j]
                 +node[verts[cn[3]]].un[0]*cell[c].b4[j]+node[verts[cn[4]]].un[0]*cell[c].b5[j]+node[verts[cn[5]]].un[0]*cell[c].b6[j]
		 +node[verts[cn[6]]].un[0]*cell[c].b7[j]+node[verts[cn[7]]].un[0]*cell[c].b8[j])/cell[c].vol;

      cell[c].dudx2[j]=(node[verts[cn[0]]].un[1]*cell[c].b1[j]+node[verts[cn[1]]].un[1]*cell[c].b2[j]+node[verts[cn[2]]].un[1]*cell[c].b3[j]
                 +node[verts[cn[3]]].un[1]*cell[c].b4[j]+node[verts[cn[4]]].un[1]*cell[c].b5[j]+node[verts[cn[5]]].un[1]*cell[c].b6[j]
                 +node[verts[cn[6]]].un[1]*cell[c].b7[j]+node[verts[cn[7]]].un[1]*cell[c].b8[j])/cell[c].vol;

      cell[c].dudx3[j]=(node[verts[cn[0]]].un[2]*cell[c].b1[j]+node[verts[cn[1]]].un[2]*cell[c].b2[j]+node[verts[cn[2]]].un[2]*cell[c].b3[j]
                 +node[verts[cn[3]]].un[2]*cell[c].b4[j]+node[verts[cn[4]]].un[2]*cell[c].b5[j]+node[verts[cn[5]]].un[2]*cell[c].b6[j]
                 +node[verts[cn[6]]].un[2]*cell[c].b7[j]+node[verts[cn[7]]].un[2]*cell[c].b8[j])/cell[c].vol;

    }

  } //end cell loop


  
}
