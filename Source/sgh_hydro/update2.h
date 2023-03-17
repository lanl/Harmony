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

void update2( 
  mesh_t  &mesh,
  mat_t &mat,
  double TIME,
  double DT,
  int NCYC,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {
/* predictor step */


  int i,j,k,l,ii,jj;

  int k1,k2,j1,j2,i1,i2;

  double E,D[10],W[10],div,mu,sum;

  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();
  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  auto cs_owned = mesh.cells_owned();

  int iorder[8],cn[8],myvert[8];
  double ave[3];

  double factor,tau10[3],tau20[3],tau30[3];

  double TOL=1.e-6;

  int numtasks,rank;
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  int kk,kk2;

  int tag=1;

  double *recv,*send;
  MPI_Status stat;    


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

  /* update nodal positions and velocities */

  for (auto v : verts ) {

    if (vflag[v]) {

      for (j=0;j<3;j++) {

	if (node[v].kfix[j]>1.e-7) node[v].fx[j]=0.0;
       
	mesh.points[v].px[j]-=node[v].dx[j];
	node[v].un[j]-=node[v].du[j];

	node[v].du[j]=DT*node[v].fx[j]/node[v].mp;
	node[v].dx[j]=DT*(node[v].un[j]+0.5*node[v].du[j]);

	mesh.points[v].px[j]+=node[v].dx[j];
	node[v].un[j]+=node[v].du[j];
      }

    }

  }

  /* calculate pdV work */
  for ( auto c : cs_owned ) {

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
    
	sum=0.0;

	//node 1

	sum-=(cell[c].f1[0]*node[verts[cn[0]]].dx[0]+cell[c].f1[1]*node[verts[cn[0]]].dx[1]+cell[c].f1[2]*node[verts[cn[0]]].dx[2])/cell[c].mass;

	//node 2

	sum-=(cell[c].f2[0]*node[verts[cn[1]]].dx[0]+cell[c].f2[1]*node[verts[cn[1]]].dx[1]+cell[c].f2[2]*node[verts[cn[1]]].dx[2])/cell[c].mass;

	//node 3

	sum-=(cell[c].f3[0]*node[verts[cn[2]]].dx[0]+cell[c].f3[1]*node[verts[cn[2]]].dx[1]+cell[c].f3[2]*node[verts[cn[2]]].dx[2])/cell[c].mass;

	//node 4

	sum-=(cell[c].f4[0]*node[verts[cn[3]]].dx[0]+cell[c].f4[1]*node[verts[cn[3]]].dx[1]+cell[c].f4[2]*node[verts[cn[3]]].dx[2])/cell[c].mass;

	//node 5

	sum-=(cell[c].f5[0]*node[verts[cn[4]]].dx[0]+cell[c].f5[1]*node[verts[cn[4]]].dx[1]+cell[c].f5[2]*node[verts[cn[4]]].dx[2])/cell[c].mass;

	//node 6

	sum-=(cell[c].f6[0]*node[verts[cn[5]]].dx[0]+cell[c].f6[1]*node[verts[cn[5]]].dx[1]+cell[c].f6[2]*node[verts[cn[5]]].dx[2])/cell[c].mass;

	//node 7

	sum-=(cell[c].f7[0]*node[verts[cn[6]]].dx[0]+cell[c].f7[1]*node[verts[cn[6]]].dx[1]+cell[c].f7[2]*node[verts[cn[6]]].dx[2])/cell[c].mass;

	//node 8

	sum-=(cell[c].f8[0]*node[verts[cn[7]]].dx[0]+cell[c].f8[1]*node[verts[cn[7]]].dx[1]+cell[c].f8[2]*node[verts[cn[7]]].dx[2])/cell[c].mass;

	cell[c].h=cell[c].h+sum;

  }


  /* update internal energy */
  for (auto c : cs_owned ) {
        
    cell[c].ec=cell[c].ec-cell[c].de+cell[c].h;

  }


  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=49;
   
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
	   //internal energy
	   send[ndat*i+index+offset]=cell[c].ec;
	   index++;	   


	   //vertex positions 
	   for (auto cr : mesh.corners_of_cell(c)) {
	     
	     auto v = mesh.corners[cr].point;

	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=mesh.points[v].px[ii];
	       index++;

	     }

	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].un[ii];
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
	   //internal energy
	   cell[c].ec=recv[ndat*i+index+offset];
	   index++;	   

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

	   //vertex positions 
	   for (int cr=0; cr<8; cr++) {
	     
	     auto v = verts[cn[cr]];

	     for (ii=0;ii<3;ii++) {
	       
               //only update vertex position for ghost nodes
	       if (vflag[v]==0) mesh.points[v].px[ii]=recv[ndat*i+index+offset];
	       index++;

	     }

	     for (ii=0;ii<3;ii++) {
	       
               //only update vertex velocity for ghost nodes
	       if (vflag[v]==0) node[v].un[ii]=recv[ndat*i+index+offset];
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
