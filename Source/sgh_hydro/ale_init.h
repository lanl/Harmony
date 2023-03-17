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
void ale_init(
  mesh_t  &mesh,
  std::vector<vertex_t> &node	      
) {


  int i,j,k,l,n;



  //save original nodal positions at begininng of Lagrangian step
  for (auto v : mesh.vertices()) { 

        //save original nodal positions at begininng of Lagrangian step
	node[v].xp[0]=node[v].x0[0];
	node[v].xp[1]=node[v].x0[1];
	node[v].xp[2]=node[v].x0[2];
  
        // zero all node-centered fluxes and gradients
	for (n=0;n<=2;n++) {
	  node[v].MUFx[n]=node[v].MVFx[n]=node[v].MWFx[n]=0.0;
	  node[v].sx_u[n]=node[v].sx_v[n]=node[v].sx_w[n]=0.0;
	}

  }




}



