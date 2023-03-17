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
void editor( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node, 
  int jedit,
  int NCYC,
  double TIME,
  double *edtimes
) {
  int i,j,k,l,id;

  FILE *out[20];

  char name[100];

  double pres,xa,ya,za,ua,va,wa,E;

  int node_id[8];

  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();

  auto cs_owned = mesh.cells_owned();

  int numtasks,rank;
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  int iorder[8],cn[8],myvert[8];
  double ave[3];

  //create ensight directory for this rank
  sprintf(name,"ls ensight_%05d",rank);

  i=system(name);

  if (i!=0) {
    sprintf(name,"mkdir ensight_%05d",rank);
    i=system(name);
    sprintf(name,"mkdir ensight_%05d/data",rank);
    i=system(name);
  }

    for ( auto c : cs ) {

      cell[c].flg=0.0;

    }


    for ( auto c : cs_owned ) {

      cell[c].flg=1.0;

    }   


  /* write geometry file */


  sprintf(name,"ensight_%05d/data/engold.%05d.geo",rank,jedit);
  out[0]=fopen(name,"w");


  fprintf(out[0],"This is the 1st description line of the EnSight Gold geometry example\n");
  fprintf(out[0],"This is the 2nd description line of the EnSight Gold geometry example\n");
  fprintf(out[0],"node id assign\n");
  fprintf(out[0],"element id assign\n");

  fprintf(out[0],"part\n");
  fprintf(out[0],"%10d\n",1);
  fprintf(out[0],"This is the description line for the part\n");
  fprintf(out[0],"coordinates\n");
  fprintf(out[0],"%10lu\n",NP);

  for (auto v : verts)
    fprintf(out[0],"%12.5e\n",mesh.points[v].px[0]);

  for (auto v : verts)
    fprintf(out[0],"%12.5e\n",mesh.points[v].px[1]);  

  for (auto v : verts)
    fprintf(out[0],"%12.5e\n",mesh.points[v].px[2]);


  fprintf(out[0],"hexa8\n");
  fprintf(out[0],"%10lu\n",NE);

  for (auto c : cs) {

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
	  node_id[iorder[i]]=v + 1;
	  i++;
	}

	fprintf(out[0],"%10lu",node_id[0]);
	fprintf(out[0],"%10lu",node_id[1]);
	fprintf(out[0],"%10lu",node_id[2]);
	fprintf(out[0],"%10lu",node_id[3]);

	fprintf(out[0],"%10lu",node_id[4]);
	fprintf(out[0],"%10lu",node_id[5]);
	fprintf(out[0],"%10lu",node_id[6]);
	fprintf(out[0],"%10lu",node_id[7]);
	fprintf(out[0],"\n");

   
  }
  fclose(out[0]);

  /* write scalar data files */
  sprintf(name,"ensight_%05d/data/engold.%05d.den",rank,jedit);
  out[0]=fopen(name,"w");
  fprintf(out[0],"Per_elem scalar values\n");
  fprintf(out[0],"part\n");
  fprintf(out[0],"%10d\n",1);
  fprintf(out[0],"hexa8\n");

  sprintf(name,"ensight_%05d/data/engold.%05d.pres",rank,jedit);
  out[1]=fopen(name,"w");
  fprintf(out[1],"Per_elem scalar values\n");
  fprintf(out[1],"part\n");
  fprintf(out[1],"%10d\n",1);
  fprintf(out[1],"hexa8\n");

  sprintf(name,"ensight_%05d/data/engold.%05d.ie",rank,jedit); 
  out[2]=fopen(name,"w");
  fprintf(out[2],"Per_elem scalar values\n");
  fprintf(out[2],"part\n");
  fprintf(out[2],"%10d\n",1);
  fprintf(out[2],"hexa8\n");

  sprintf(name,"ensight_%05d/data/engold.%05d.flg",rank,jedit); 
  out[3]=fopen(name,"w");
  fprintf(out[3],"Per_elem scalar values\n");
  fprintf(out[3],"part\n");
  fprintf(out[3],"%10d\n",1);
  fprintf(out[3],"hexa8\n");

  sprintf(name,"ensight_%05d/data/engold.%05d.vol",rank,jedit); 
  out[4]=fopen(name,"w");
  fprintf(out[4],"Per_elem scalar values\n");
  fprintf(out[4],"part\n");
  fprintf(out[4],"%10d\n",1);
  fprintf(out[4],"hexa8\n");

  
  for (auto c : cs) {

    pres=cell[c].pc;

    if (pres>0. && pres<1.e-32) pres=0.0;
    if (pres<0. && pres>-1.e-32) pres=0.0;

    fprintf(out[0],"%12.5e\n",cell[c].dc);
    fprintf(out[1],"%12.5e\n",pres);
    fprintf(out[2],"%12.5e\n",cell[c].ec);
    fprintf(out[3],"%12.5e\n",cell[c].flg);
    fprintf(out[4],"%12.5e\n",cell[c].vol);
    
  }
  
  for (i=0;i<5;i++)
    fclose(out[i]);

  /* write vector data files */  
  sprintf(name,"ensight_%05d/data/engold.%05d.vel",rank,jedit);

  out[0]=fopen(name,"w");
  fprintf(out[0],"Per_node vector values\n");
  fprintf(out[0],"part\n");
  fprintf(out[0],"%10d\n",1);
  fprintf(out[0],"coordinates\n");

  for (auto v : verts) {
    E=node[v].un[0];
    if (E<0.0) E=-1.0*E;
    if (E<1.e-32)
      fprintf(out[0],"%12.5e\n",0.0);
    else
      fprintf(out[0],"%12.5e\n",node[v].un[0]);
  }

  for (auto v : verts) {
    E=node[v].un[1];
    if (E<0.0) E=-1.0*E;
    if (E<1.e-32)
      fprintf(out[0],"%12.5e\n",0.0);
    else
      fprintf(out[0],"%12.5e\n",node[v].un[1]);
  }

  for (auto v : verts) {
    E=node[v].un[2];
    if (E<0.0) E=-1.0*E;
    if (E<1.e-32)
      fprintf(out[0],"%12.5e\n",0.0);
    else
      fprintf(out[0],"%12.5e\n",node[v].un[2]);
  }  

  fclose(out[0]);



  /* write vector boundary specifications */  
  sprintf(name,"ensight_%05d/data/engold.%05d.kfix",rank,jedit);

  out[0]=fopen(name,"w");
  fprintf(out[0],"Per_node vector values\n");
  fprintf(out[0],"part\n");
  fprintf(out[0],"%10d\n",1);
  fprintf(out[0],"coordinates\n");

  for (auto v : verts) {
    auto & kfx = node[v].kfix;
    fprintf(out[0],"%12.5e\n",kfx[0]);
  }

  for (auto v : verts)  {
    auto & kfx = node[v].kfix;
    fprintf(out[0],"%12.5e\n",kfx[1]); 
  } 

  for (auto v : verts)  {
    auto & kfx = node[v].kfix;
    fprintf(out[0],"%12.5e\n",kfx[2]);
  }

  fclose(out[0]);

  
  sprintf(name,"ensight_%05d/data/time.dat",rank);
  out[0]=fopen(name,"w");

  fprintf(out[0],"%d %d %12.5e\n",jedit+1,NCYC,TIME);

  fclose(out[0]);


}
