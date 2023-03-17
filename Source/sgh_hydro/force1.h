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
void force1( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon,
  double DT,
  double C1,
  double C2,
  double C3
) {

  int i,j,k,l;

  int k1,k2,j1,j2,i1,i2;

  double E, tau[10], D[10];

  double factor;
  
  double areab[25],ac[25],force[25],muc[9],sum[5],ave[4];

  double q,DVC,length;

  double ua,va,wa;
  double xa,ya,za;
  double dx,dy,dz;


  
  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();
  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  auto cs_owned = mesh.cells_owned();

  int rank,numtasks,kk,kk2,ii;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  MPI_Status stat;

  double *recv,*send;
  

  int iorder[8],cn[8],myvert[8];

  int *vflag;

  double corner_force[3];

  double phi;



  // loop over cells to find corner forces
  for (auto c : cs_owned ) {

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

	// calculate the Cauchy stress tensor 
	// first, copy the stress deviators
	tau(1,1)=cell[c].stress1[0];
	tau(1,2)=cell[c].stress1[1];
	tau(1,3)=cell[c].stress1[2];
	tau(2,1)=cell[c].stress2[0];
	tau(2,2)=cell[c].stress2[1];
	tau(2,3)=cell[c].stress2[2];
	tau(3,1)=cell[c].stress3[0];
	tau(3,2)=cell[c].stress3[1];
	tau(3,3)=cell[c].stress3[2];

	//find artifical viscosity
	DVC=cell[c].dudx1[0]+cell[c].dudx2[1]+cell[c].dudx3[2];
	length=pow(cell[c].vol,1./3.);

	if (DVC<0.0) 
	  q=cell[c].dc*(C1*length*length*DVC*DVC-C2*length*cell[c].sound*DVC);
	else
	  q=0.;	
	
	// then add pressure and artificial viscosity
	tau(1,1)-=(cell[c].pc+q);
	tau(2,2)-=(cell[c].pc+q);
	tau(3,3)-=(cell[c].pc+q);

	// finally reverse the sign
	tau(1,1)*=-1.0;
	tau(1,2)*=-1.0;
	tau(1,3)*=-1.0;
	tau(2,1)*=-1.0;
	tau(2,2)*=-1.0;
	tau(2,3)*=-1.0;
	tau(3,1)*=-1.0;
	tau(3,2)*=-1.0;
	tau(3,3)*=-1.0;

	// zero pDV work term
	cell[c].h=0.0;
    
    
	// surface areas associated with each vertex
	areab(1,1)=cell[c].b1[0];
	areab(2,1)=cell[c].b1[1];
	areab(3,1)=cell[c].b1[2];

	areab(1,2)=cell[c].b2[0];
	areab(2,2)=cell[c].b2[1];
	areab(3,2)=cell[c].b2[2];

	areab(1,3)=cell[c].b3[0];
	areab(2,3)=cell[c].b3[1];
	areab(3,3)=cell[c].b3[2];

	areab(1,4)=cell[c].b4[0];
	areab(2,4)=cell[c].b4[1];
	areab(3,4)=cell[c].b4[2];

	areab(1,5)=cell[c].b5[0];
	areab(2,5)=cell[c].b5[1];
	areab(3,5)=cell[c].b5[2];

	areab(1,6)=cell[c].b6[0];
	areab(2,6)=cell[c].b6[1];
	areab(3,6)=cell[c].b6[2];

	areab(1,7)=cell[c].b7[0];
	areab(2,7)=cell[c].b7[1];
	areab(3,7)=cell[c].b7[2];

	areab(1,8)=cell[c].b8[0];
	areab(2,8)=cell[c].b8[1];
	areab(3,8)=cell[c].b8[2];


	//find the corner forces
	//node 1
	cell[c].f1[0]=tau(1,1)*areab(1,1)+tau(1,2)*areab(2,1)+tau(1,3)*areab(3,1); 
	cell[c].f1[1]=tau(2,1)*areab(1,1)+tau(2,2)*areab(2,1)+tau(2,3)*areab(3,1); 
	cell[c].f1[2]=tau(3,1)*areab(1,1)+tau(3,2)*areab(2,1)+tau(3,3)*areab(3,1);

	//node 2
	cell[c].f2[0]=tau(1,1)*areab(1,2)+tau(1,2)*areab(2,2)+tau(1,3)*areab(3,2);
	cell[c].f2[1]=tau(2,1)*areab(1,2)+tau(2,2)*areab(2,2)+tau(2,3)*areab(3,2);
	cell[c].f2[2]=tau(3,1)*areab(1,2)+tau(3,2)*areab(2,2)+tau(3,3)*areab(3,2);

	//node 3
	cell[c].f3[0]=tau(1,1)*areab(1,3)+tau(1,2)*areab(2,3)+tau(1,3)*areab(3,3);
	cell[c].f3[1]=tau(2,1)*areab(1,3)+tau(2,2)*areab(2,3)+tau(2,3)*areab(3,3);
	cell[c].f3[2]=tau(3,1)*areab(1,3)+tau(3,2)*areab(2,3)+tau(3,3)*areab(3,3);

	//node 4
	cell[c].f4[0]=tau(1,1)*areab(1,4)+tau(1,2)*areab(2,4)+tau(1,3)*areab(3,4);
	cell[c].f4[1]=tau(2,1)*areab(1,4)+tau(2,2)*areab(2,4)+tau(2,3)*areab(3,4);
	cell[c].f4[2]=tau(3,1)*areab(1,4)+tau(3,2)*areab(2,4)+tau(3,3)*areab(3,4);

	//node 5
	cell[c].f5[0]=tau(1,1)*areab(1,5)+tau(1,2)*areab(2,5)+tau(1,3)*areab(3,5);
	cell[c].f5[1]=tau(2,1)*areab(1,5)+tau(2,2)*areab(2,5)+tau(2,3)*areab(3,5);
	cell[c].f5[2]=tau(3,1)*areab(1,5)+tau(3,2)*areab(2,5)+tau(3,3)*areab(3,5);

	//node 6
	cell[c].f6[0]=tau(1,1)*areab(1,6)+tau(1,2)*areab(2,6)+tau(1,3)*areab(3,6);
	cell[c].f6[1]=tau(2,1)*areab(1,6)+tau(2,2)*areab(2,6)+tau(2,3)*areab(3,6);
	cell[c].f6[2]=tau(3,1)*areab(1,6)+tau(3,2)*areab(2,6)+tau(3,3)*areab(3,6);

	//node 7
	cell[c].f7[0]=tau(1,1)*areab(1,7)+tau(1,2)*areab(2,7)+tau(1,3)*areab(3,7);
	cell[c].f7[1]=tau(2,1)*areab(1,7)+tau(2,2)*areab(2,7)+tau(2,3)*areab(3,7);
	cell[c].f7[2]=tau(3,1)*areab(1,7)+tau(3,2)*areab(2,7)+tau(3,3)*areab(3,7);

	//node 8
	cell[c].f8[0]=tau(1,1)*areab(1,8)+tau(1,2)*areab(2,8)+tau(1,3)*areab(3,8);
	cell[c].f8[1]=tau(2,1)*areab(1,8)+tau(2,2)*areab(2,8)+tau(2,3)*areab(3,8);
	cell[c].f8[2]=tau(3,1)*areab(1,8)+tau(3,2)*areab(2,8)+tau(3,3)*areab(3,8);

  

	
  } //end cell loop


  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=24;
   
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

	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f1[ii];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f2[ii];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f3[ii];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f4[ii];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f5[ii];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f6[ii];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f7[ii];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     send[ndat*i+index+offset]=cell[c].f8[ii];
	     index++;
	     
	   }	   


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

	   //loop over all corners in the cell
	   int index=0;

	   for (ii=0;ii<3;ii++) {

	     cell[c].f1[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     cell[c].f2[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     cell[c].f3[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     cell[c].f4[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     cell[c].f5[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     cell[c].f6[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }
	   for (ii=0;ii<3;ii++) {

	     cell[c].f7[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }	   
	   for (ii=0;ii<3;ii++) {

	     cell[c].f8[ii]=recv[ndat*i+index+offset];
	     index++;
	     
	   }	   	   
  

	 }//end cell loop

       } //end if recv_num[kk]>0
	 
       
       kk2++;

     } //end if rank!=j

   } //done loop over ranks


   delete[] send;
   delete[] recv;


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
    
   // loop over all vertices and accumulate the corner forces
  for (auto v : verts )  {
    
    i=v;
    if (vflag[i]) {
      
    node[v].fx[0]=node[v].fx[1]=node[v].fx[2]=0.0;

    //   if (rank==1) printf("%d %f %f %f\n",j,mesh.points[v].px[0],mesh.points[v].px[1],mesh.points[v].px[2]);

    // loop over cells associated with this vertex
    for (auto c : mesh.cells_of_point(v) ) {



        myvert[0]=cell[c].iord1;
        myvert[1]=cell[c].iord2;
        myvert[2]=cell[c].iord3;
        myvert[3]=cell[c].iord4;
        myvert[4]=cell[c].iord5;
        myvert[5]=cell[c].iord6;
        myvert[6]=cell[c].iord7;
        myvert[7]=cell[c].iord8;

	j=0;
	for ( auto v2 : mesh.vertices_of_cell(c) ) {

	  for (i=0;i<8;i++) 
            if (node[v2].gid==myvert[i]) iorder[j]=i;

	  j++;
	}	

	

      //find the corner that is associted with v
      i=l=0;
      for ( auto v2 : mesh.vertices_of_cell(c) ) {
	k=node[v2].gid;
	j=node[v].gid;
	if (j==k) l=i;
        i++;
      }
      
      // if (rank==0) printf("%d %d\n",j,l);
      
      l=iorder[l];
	corner_force[0] = cell[c].f1[0];
	corner_force[1] = cell[c].f1[1];
	corner_force[2] = cell[c].f1[2];
      
	if (l==1) {
	  corner_force[0] = cell[c].f2[0];
	  corner_force[1] = cell[c].f2[1];
	  corner_force[2] = cell[c].f2[2];	
	}  
	else if (l==2) {
	  corner_force[0] = cell[c].f3[0];
	  corner_force[1] = cell[c].f3[1];
	  corner_force[2] = cell[c].f3[2];	
	}
	else if (l==3) {
	  corner_force[0] = cell[c].f4[0];
	  corner_force[1] = cell[c].f4[1];
	  corner_force[2] = cell[c].f4[2];	
	}
	else if (l==4) {
	  corner_force[0] = cell[c].f5[0];
	  corner_force[1] = cell[c].f5[1];
	  corner_force[2] = cell[c].f5[2];	
	}
	else if (l==5) {
	  corner_force[0] = cell[c].f6[0];
	  corner_force[1] = cell[c].f6[1];
	  corner_force[2] = cell[c].f6[2];	
	}
	else if (l==6) {
	  corner_force[0] = cell[c].f7[0];
	  corner_force[1] = cell[c].f7[1];
	  corner_force[2] = cell[c].f7[2];	
	}
	else if (l==7) {
	  corner_force[0] = cell[c].f8[0];
	  corner_force[1] = cell[c].f8[1];
	  corner_force[2] = cell[c].f8[2];	
	}
      
	node[v].fx[0]+=corner_force[0];
	node[v].fx[1]+=corner_force[1];    
	node[v].fx[2]+=corner_force[2];
      
      }  //end cell loop

    } //end vflag test
	 
  } // end vertex loop

  delete[] vflag;


}
