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

void read_connect(
     mesh_connect_t &mcon
){

  int i,j,k,l,id,m,n,p,kk,flg2;

  FILE *in;
  
  int local;
  int remote;

  char name[100];
  int rank,count;

  int num_send,num_recv;


  //get the rank and the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &m);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //read the connectivity file for this rank to find the number of cells sent and received by this rank
  sprintf(name,"mesh/ensight_%05d/connectivity.dat",rank);
  in=fopen(name,"r");
  //read the send array from the connectivity file
  count=0;
  for (j=0;j<m;j++) { //loop over all ranks
    fscanf(in,"%d",&k); //the number of cells that this rank sends to rank j
    
    count+=k;

    //read the local and remote ids of the cells that this rank sends to rank j
    for (i=0;i<k;i++) {
      fscanf(in,"%d%d",&local, &remote);
    }

  }
  num_send=count;
  count=0;
  for (j=0;j<m;j++) { //loop over all ranks
    fscanf(in,"%d",&k); //the number of cells that rank j sends to this rank
    
    count+=k;
    
    //read the local and remote ids of the  cells that rank j sends to this rank
    for (i=0;i<k;i++) {
      fscanf(in,"%d%d",&local, &remote);
    }

  }
  num_recv=count;
  
  fclose(in); //done reading connectivity file

  //size the mesh.send and mesh.recv arrays for this rank
  mcon.send_local.reserve(num_send);
  mcon.send_rank.reserve(num_send);  
  mcon.send_cell.reserve(num_send);

  mcon.send_local.resize(0);
  mcon.send_rank.resize(0);  
  mcon.send_cell.resize(0);

  mcon.send_num.reserve(m-1);
  mcon.send_num.resize(0);  
  
  mcon.recv_local.reserve(num_recv);
  mcon.recv_rank.reserve(num_recv);  
  mcon.recv_cell.reserve(num_recv);

  mcon.recv_local.resize(0);
  mcon.recv_rank.resize(0);  
  mcon.recv_cell.resize(0);

  mcon.recv_num.reserve(m-1);
  mcon.recv_num.resize(0);  
  
  //read the connectivity file a second time for this rank
  sprintf(name,"mesh/ensight_%05d/connectivity.dat",rank);
  in=fopen(name,"r");
  //read the send array from the connectivity file
  for (j=0;j<m;j++) { //loop over all ranks
    fscanf(in,"%d",&k); //the number of cells that this rank sends to rank j
    
    if (j!=rank) mcon.send_num.push_back(k);

    //read the local and remote ids of the cells that this rank sends to rank j
    for (i=0;i<k;i++) {
      fscanf(in,"%d%d",&local, &remote);
      mcon.send_local.push_back(local);
      mcon.send_cell.push_back(remote);
      mcon.send_rank.push_back(j);
    }

  }
  //read the recv array from the connectivity file
  for (j=0;j<m;j++) { //loop over all ranks
    fscanf(in,"%d",&k); //the number of cells that this rank sends to rank j

    if (j!=rank) mcon.recv_num.push_back(k);

    //read the local and remote ids of the cells that this rank sends to rank j
    for (i=0;i<k;i++) {
      fscanf(in,"%d%d",&local, &remote);
      mcon.recv_local.push_back(local);
      mcon.recv_cell.push_back(remote);
      mcon.recv_rank.push_back(j);
    }

  }
  
  fclose(in); //done reading connectivity file
  /*
  if (rank==0) {

    printf("send array\n");
    for (i=0;i<num_send;i++) {

      printf("%d %d\n",mcon.send_local[i],mcon.send_cell[i]);

    }

    printf("receive array\n");
    for (i=0;i<num_recv;i++) {

      printf("%d %d\n",mcon.recv_local[i],mcon.recv_cell[i]);

    }
  }
 
  */
}

