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
void ale_update( 
  mesh_t  &mesh,
  mat_t &mat,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon
) {

  int i,j,k,p,l;

  int imin,imax,jmin,jmax,kmin,kmax;

  int n1,i0,i1,i2,i3,i4,i5,i6;

  double vol,*volp,mflux,IE,EN,mass,mp1;
  double dx,dy,dz;
  double x,y,z;
  double sum,sum1,factor,Vm;

  double x1,x2,x3,x4,x5,x6,x7,x8;
  double y1,y2,y3,y4,y5,y6,y7,y8;
  double z1,z2,z3,z4,z5,z6,z7,z8;


  //vertices
  auto verts = mesh.vertices();
  auto vp = verts[0];

  //cells
  auto cs = mesh.cells();
  auto c1 = cs[0];

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

  /* update density and ie in all cells */
  for ( auto c : mesh.cells_owned() ) {


    /* indicies of vertex points */
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


    i=c;
    /* indicies of face neighbors */
    for (j=0;j<6;j++)
      indx[j]=-1; //default indx of -1 indicates that the cell is on boundary

    //find +x face
    myvert[0]=cn[1];
    myvert[1]=cn[2];
    myvert[2]=cn[5];
    myvert[3]=cn[6];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++)
	  if (myvert[j]==v1) counter++;
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[0]=j;
    }

    //find -x face
    myvert[0]=cn[0];
    myvert[1]=cn[3];
    myvert[2]=cn[4];
    myvert[3]=cn[7];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++)
	  if (myvert[j]==v1) counter++;
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[1]=j;
    }

    //find +y face
    myvert[0]=cn[3];
    myvert[1]=cn[2];
    myvert[2]=cn[6];
    myvert[3]=cn[7];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++)
	  if (myvert[j]==v1) counter++;
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[2]=j;
    }

    //find -y face
    myvert[0]=cn[0];
    myvert[1]=cn[1];
    myvert[2]=cn[4];
    myvert[3]=cn[5];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++)
	  if (myvert[j]==v1) counter++;
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[3]=j;
    }

    //find +z face
    myvert[0]=cn[4];
    myvert[1]=cn[5];
    myvert[2]=cn[6];
    myvert[3]=cn[7];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++)
	  if (myvert[j]==v1) counter++;
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[4]=j;
    }

    //find -z face
    myvert[0]=cn[0];
    myvert[1]=cn[1];
    myvert[2]=cn[2];
    myvert[3]=cn[3];

    for (auto f : mesh.faces_of_cell(c)) {
      counter=0;
      for (auto v1 : mesh.points_of_face(f))
	for (j=0;j<4;j++) {
	  if (myvert[j]==v1) counter++;
	}
      if (counter==4) myface=f;
    }

    for (auto c1 : mesh.cells_of_face(myface)) {
      j=c1;
      if (j!=i) indx[5]=j;
    }


    
    IE=cell[c].ec*cell[c].mass;
    EN=cell[c].enc*cell[c].mass;
    mass=cell[c].mass;

    //add flux contribution from +x face
    if (indx[0]!=-1) {
      mass-=cell[c].MFx[0];
      IE-=cell[c].EFx[0];
      EN-=cell[c].TFx[0];
    }

    //add flux contribution from -x face
    if (indx[1]!=-1) {

      c1=cs[indx[1]];
      mass+=cell[c1].MFx[0];
      IE+=cell[c1].EFx[0];
      EN+=cell[c1].TFx[0];
    }

    //add flux contribution from +y face
    if (indx[2]!=-1) {
      mass-=cell[c].MFx[1];
      IE-=cell[c].EFx[1];
      EN-=cell[c].TFx[1];
    }

    //add flux contribution from -y face
    if (indx[3]!=-1) {

      c1=cs[indx[3]];
      mass+=cell[c1].MFx[1];
      IE+=cell[c1].EFx[1];
      EN+=cell[c1].TFx[1];
    }

    //add flux contribution from +z face
    if (indx[4]!=-1) {
	
      mass-=cell[c].MFx[2];
      IE-=cell[c].EFx[2];
      EN-=cell[c].TFx[2];
    }

    //add flux contribution from -z face
    if (indx[5]!=-1) {

      c1=cs[indx[5]];
      mass+=cell[c1].MFx[2];
      IE+=cell[c1].EFx[2];
      EN+=cell[c1].TFx[2];
    }

    //update mass, internal energy and total energy 
    cell[c].mass=mass;
    cell[c].ec=IE/(mass+1.e-15);
    cell[c].enc=EN/(mass+1.e-15);

    // find new cell volume
    x1=mesh.points[cn[0]].px[0];
    x2=mesh.points[cn[1]].px[0];
    x3=mesh.points[cn[2]].px[0];
    x4=mesh.points[cn[3]].px[0];
    x5=mesh.points[cn[4]].px[0];
    x6=mesh.points[cn[5]].px[0];
    x7=mesh.points[cn[6]].px[0];
    x8=mesh.points[cn[7]].px[0];

    y1=mesh.points[cn[0]].px[1];
    y2=mesh.points[cn[1]].px[1];
    y3=mesh.points[cn[2]].px[1];
    y4=mesh.points[cn[3]].px[1];
    y5=mesh.points[cn[4]].px[1];
    y6=mesh.points[cn[5]].px[1];
    y7=mesh.points[cn[6]].px[1];
    y8=mesh.points[cn[7]].px[1];

    z1=mesh.points[cn[0]].px[2];
    z2=mesh.points[cn[1]].px[2];
    z3=mesh.points[cn[2]].px[2];
    z4=mesh.points[cn[3]].px[2];
    z5=mesh.points[cn[4]].px[2];
    z6=mesh.points[cn[5]].px[2];
    z7=mesh.points[cn[6]].px[2];
    z8=mesh.points[cn[7]].px[2]; 

    /* new cell volume */
    vol= (x2*(y4*(-z1 + z3) + y5*(z1 - z6) + y1*(z3 + z4 - z5 - z6) + y7*(-z3 + z6) + y6*(z1 - z3 + z5 - z7) + y3*(-z1 - z4 + z6 + z7)) +
               x8*(y1*(-z4 + z5) + y7*(z3 + z4 - z5 - z6) + y3*(z4 - z7) + y4*(z1 - z3 + z5 - z7) + y6*(-z5 + z7) + y5*(-z1 - z4 + z6 + z7)) +
               x4*(y2*(z1 - z3) + y8*(-z1 + z3 - z5 + z7) + y7*(z3 - z8) + y3*(z1 + z2 - z7 - z8) + y5*(-z1 + z8) + y1*(-z2 - z3 + z5 + z8)) + 
               x6*(y1*(z2 - z5) + y8*(z5 - z7) + y3*(-z2 + z7) + y2*(-z1 + z3 - z5 + z7) + y5*(z1 + z2 - z7 - z8) + y7*(-z2 - z3 + z5 + z8)) + 
               x7*(y2*(z3 - z6) + y8*(-z3 - z4 + z5 + z6) + y6*(z2 + z3 - z5 - z8) + y5*(z6 - z8) + y4*(-z3 + z8) + y3*(-z2 + z4 - z6 + z8)) + 
               x1*(y3*(z2 - z4) + y8*(z4 - z5) + y6*(-z2 + z5) + y2*(-z3 - z4 + z5 + z6) + y4*(z2 + z3 - z5 - z8) + y5*(-z2 + z4 - z6 + z8)) + 
               x3*(y1*(-z2 + z4) + y6*(z2 - z7) + y2*(z1 + z4 - z6 - z7) + y8*(-z4 + z7) + y7*(z2 - z4 + z6 - z8) + y4*(-z1 - z2 + z7 + z8)) + 
	      x5*(y2*(-z1 + z6) + y8*(z1 + z4 - z6 - z7) + y4*(z1 - z8) + y1*(z2 - z4 + z6 - z8) + y7*(-z6 + z8) + y6*(-z1 - z2 + z7 + z8)))/12.;    

  } //end cell loop



  //MPI communication with other ranks
  //construct array to hold data to send to other ranks
  
  int ndat=3;
   
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

	   send[ndat*i+index+offset]=cell[c].mass;
	   index++;
	   send[ndat*i+index+offset]=cell[c].ec;
	   index++;
	   send[ndat*i+index+offset]=cell[c].enc;
	   index++;
	   
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

	   cell[c].mass=recv[ndat*i+index+offset];
	   index++;
	   cell[c].ec=recv[ndat*i+index+offset];
	   index++;
	   cell[c].enc=recv[ndat*i+index+offset];
	   index++;

	 }//end cell loop

       } //end if recv_num[kk]>0
	 
       
       kk2++;

     } //end if rank!=j

   } //done loop over ranks


   delete[] send;
   delete[] recv;
  
}


