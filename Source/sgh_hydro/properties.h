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

void properties( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell 
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
  
  int iorder[8],cn[8],myvert[8];
  double ave[3];

  double factor,tau10[3],tau20[3],tau30[3];

  const double TOL=1.e-6;

  /* calculate the cell volume, density and forces at the vertices */
  for ( auto c : cs ) {

    /* stress tensor components set to zero */
    cell[c].stress1[0]=cell[c].stress1[1]=cell[c].stress1[2]=0.0;
    cell[c].stress2[0]=cell[c].stress2[1]=cell[c].stress2[2]=0.0;
    cell[c].stress3[0]=cell[c].stress3[1]=cell[c].stress3[2]=0.0;
    
    /* cell properties */
    cell[c].dc=cell[c].mass/cell[c].vol;
    cell[c].pc=gas1(mat,cell[c].dc,cell[c].ec);
    cell[c].sound=gas2(mat,cell[c].dc,cell[c].ec);
    
  } //end cell loop


  
}
