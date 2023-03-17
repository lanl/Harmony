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
void ale_update2( 
  mesh_t  &mesh,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {


  int i,j,k,p,l,ii,jj;

  int imin,imax,jmin,jmax,kmin,kmax;

  int n1,i0,i1,i2,i3,i4,i5,i6;

  double vol,volp,mflux,IE,EN,mass,mp1;
  double dx,dy,dz;
  double x,y,z;
  double sum,sum1;

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;


  //vertices
  auto verts = mesh.vertices();
  auto vp = verts[0];
  auto NP = verts.size();
  
  //cells
  auto cs = mesh.cells();
  auto c1 = cs[0];

  auto cs_owned = mesh.cells_owned();

  auto faces = mesh.list_faces();
  auto myface = faces[0];

  int iorder[8],cn[8],myvert[8],counter,indx[27];

  int rank,numtasks,kk,kk2,k2;

  int tag=1;
   
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
 
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

  /* update u,v,w velocities */
  for ( auto v : verts ) {

    if (vflag[v]) {

      // temporary storage of the new vertex mass in mp1
      mp1=0.0;
      //gather corner masses
      for ( auto c : mesh.cells_of_point(v) ) 
	mp1+=0.125*cell[c].mass;
	
      // indicies of neighboring nodes 
      //begin construction of the indx array containing vertex neighbors
      i=v;
      for (j=0;j<27;j++)
	indx[j]=-1; //default setting of -1 indicates boundary vertex

      indx[13]=v; //central vertex on the 27 point stencil

      //loop over all cells containing vertex v
      for ( auto c : mesh.cells_of_point(v) ) {

	//order the vertices in the cell c
	myvert[0]=cell[c].iord1;
	myvert[1]=cell[c].iord2;
	myvert[2]=cell[c].iord3;
	myvert[3]=cell[c].iord4;
	myvert[4]=cell[c].iord5;
	myvert[5]=cell[c].iord6;
	myvert[6]=cell[c].iord7;
	myvert[7]=cell[c].iord8;

	j=0;
	for ( auto v1 : mesh.vertices_of_cell(c) ) {
	  for (k=0;k<8;k++)
	    if (node[v1].gid==myvert[k]) iorder[j]=k;

	  j++;
	}

    
	k=0;
	for ( auto v1 : mesh.vertices_of_cell(c) ) {
	  cn[iorder[k]]=v1;
	  k++;
	}

	//find which vertex corresponds to v
	for (k=0,l=0;k<8;k++) {
	  if (cn[k]==i) l=k;
	}


	if (l==0) {

	  indx[13]=cn[0];
	  indx[14]=cn[1];
	  indx[17]=cn[2];
	  indx[16]=cn[3];
	  indx[22]=cn[4];
	  indx[23]=cn[5];
	  indx[26]=cn[6];
	  indx[25]=cn[7];

	}
	else if (l==1) {

	  indx[12]=cn[0];
	  indx[13]=cn[1];
	  indx[16]=cn[2];
	  indx[15]=cn[3];
	  indx[21]=cn[4];
	  indx[22]=cn[5];
	  indx[25]=cn[6];
	  indx[24]=cn[7];

	}

	else if (l==2) {

	  indx[9]=cn[0];
	  indx[10]=cn[1];
	  indx[13]=cn[2];
	  indx[12]=cn[3];
	  indx[18]=cn[4];
	  indx[19]=cn[5];
	  indx[22]=cn[6];
	  indx[21]=cn[7];


	}
	else if (l==3) {

	  indx[10]=cn[0];
	  indx[11]=cn[1];
	  indx[14]=cn[2];
	  indx[13]=cn[3];
	  indx[19]=cn[4];
	  indx[20]=cn[5];
	  indx[23]=cn[6];
	  indx[22]=cn[7];


	}

	else if (l==4) {

	  indx[4]=cn[0];
	  indx[5]=cn[1];
	  indx[8]=cn[2];
	  indx[7]=cn[3];
	  indx[13]=cn[4];
	  indx[14]=cn[5];
	  indx[17]=cn[6];
	  indx[16]=cn[7];


	}
	else if (l==5) {

	  indx[3]=cn[0];
	  indx[4]=cn[1];
	  indx[7]=cn[2];
	  indx[6]=cn[3];
	  indx[12]=cn[4];
	  indx[13]=cn[5];
	  indx[16]=cn[6];
	  indx[15]=cn[7];


	}
	else if (l==6) {

	  indx[0]=cn[0];
	  indx[1]=cn[1];
	  indx[4]=cn[2];
	  indx[3]=cn[3];
	  indx[9]=cn[4];
	  indx[10]=cn[5];
	  indx[13]=cn[6];
	  indx[12]=cn[7];


	}
	else if (l==7) {

	  indx[1]=cn[0];
	  indx[2]=cn[1];
	  indx[5]=cn[2];
	  indx[4]=cn[3];
	  indx[10]=cn[4];
	  indx[11]=cn[5];
	  indx[14]=cn[6];
	  indx[13]=cn[7];


	}

      } //finished looping over all cells containing vertex v
      //indx array is now built

      //begin velocity update
      node[v].un[0]=node[v].un[0]*node[v].mp;
      node[v].un[1]=node[v].un[1]*node[v].mp;
      node[v].un[2]=node[v].un[2]*node[v].mp;

      //add flux contribution from +x face
      if (indx[14]!=-1) {

	node[v].un[0]-=node[v].MUFx[0];
	node[v].un[1]-=node[v].MVFx[0];
	node[v].un[2]-=node[v].MWFx[0];

      }

      //add flux contribution from -x face
      if (indx[12]!=-1) {

	vp=verts[indx[12]];
	node[v].un[0]+=node[vp].MUFx[0];
	node[v].un[1]+=node[vp].MVFx[0];
	node[v].un[2]+=node[vp].MWFx[0];

      }

      //add flux contribution from +y face
      if (indx[16]!=-1) {

	node[v].un[0]-=node[v].MUFx[1];
	node[v].un[1]-=node[v].MVFx[1];
	node[v].un[2]-=node[v].MWFx[1];

      }

      //add flux contribution from -y face
      if (indx[10]!=-1) {

	vp=verts[indx[10]];
	node[v].un[0]+=node[vp].MUFx[1];
	node[v].un[1]+=node[vp].MVFx[1];
	node[v].un[2]+=node[vp].MWFx[1];

      }

      //add flux contribution from +z face
      if (indx[22]!=-1) {

	node[v].un[0]-=node[v].MUFx[2];
	node[v].un[1]-=node[v].MVFx[2];
	node[v].un[2]-=node[v].MWFx[2];

      }

      //add flux contribution from -z face
      if (indx[4]!=-1) {

	vp=verts[indx[4]];
	node[v].un[0]+=node[vp].MUFx[2];
	node[v].un[1]+=node[vp].MVFx[2];
	node[v].un[2]+=node[vp].MWFx[2];

      }

      //update vertex velocity
      node[v].un[0]/=mp1;
      node[v].un[1]/=mp1;
      node[v].un[2]/=mp1;

      //update vertex mass
      node[v].mp=mp1;

      //apply velocity boundary conditions
      for (j=0;j<3;j++) {

	if (node[v].kfix[j]>1.e-7) node[v].un[j]=0.0;
      }

    } // if vflag equals 1


  } //end vertex loop


  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=32;
   
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

	   //loop over all corners in the cell
	   for (auto cr : mesh.corners_of_cell(c)) {
	     
	     auto v = mesh.corners[cr].point;

	     //nodal velocities
	     for (ii=0;ii<3;ii++) {
	      
	       send[ndat*i+index+offset]=node[v].un[ii];
	       index++;

	     }
	     //nodal mass
	     send[ndat*i+index+offset]=node[v].mp;
	     index++;     
	     

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

	   //loop over corners 
	   for (int cr=0; cr<8; cr++) {
	     
	     auto v = verts[cn[cr]];

	     //nodal velocities
	     for (ii=0;ii<3;ii++) {
	       
               //only update nodal velocities for ghost nodes
	       if (vflag[v]==0) node[v].un[ii]=recv[ndat*i+index+offset];
	       index++;

	     }
	     //nodal mass
	     //only update ghost nodes
	     if (vflag[v]==0) node[v].mp=recv[ndat*i+index+offset];
	     index++;

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
