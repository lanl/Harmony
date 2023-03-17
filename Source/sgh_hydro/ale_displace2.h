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
void ale_displace2( 
  mesh_t  &mesh,
  std::vector<vertex_t> &node
) {

  /* set xp,yp,zp positions equal to positions after Lagrange step */
  for (auto v : mesh.vertices()) { 
 
    node[v].xp[0]=mesh.points[v].px[0];
    node[v].xp[1]=mesh.points[v].px[1];
    node[v].xp[2]=mesh.points[v].px[2];

  }

}
