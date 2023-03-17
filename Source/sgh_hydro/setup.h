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
class trilink_cell_t {
  
public :

  double Vc;
  double Mc;
  double dc;
  double ec;
  double vel[2];
  
  trilink_cell_t();

};  
  
trilink_cell_t::trilink_cell_t()
{

  Vc=0.0;
  Mc=0.0;
  dc=0.0;
  ec=0.0;

  vel[0]=vel[1]=0.0;

}



void setup( 
  mesh_t  &mesh,
  mat_t &mat,
  global_data_t &gd,
  fill_data_t &fd,
  bound_data_t &bd,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node,
  mesh_connect_t &mcon  
) {


  int i,j,k,l,ibkg,ii,jj,kk;
  int i1,i2,j1,j2,k1,k2,JT;

  double dx,dy,dz;
  double x[8],y[8],z[8];
  double x1,x2;
  double y1,y2;
  double z1,z2;

  double rho,E;

  double theta1, theta2, rad2,theta, dtheta;
  double phi1, phi2, phi, dphi;
  double avex,avey,avez;
 
  FILE *in;

  double *th1,*r1,*yp1,rad,Pi=2.*acos(0.);
  int num;

  int ix,iy,iz,cell_id;

  //vertices
  auto verts = mesh.vertices();
  auto NP = verts.size();

  //cells
  auto cs = mesh.cells();
  auto NE = cs.size();

  auto cs_owned = mesh.cells_owned();
          
  int iorder[8],cn[8],myvert[8];
  double ave[3],ave2[3],dx1[3],dx2[3],cross[3];

  int *pn,kkpll,kkzln;
  double *px, *py, *pz, *ap,xa;

  char ch;

  double maxv[3],minv[3];

  int jedit,mid;
  char name[100];
  FILE *out[100];
  double TOL;

  TOL=1.e-6;

  int numtasks,rank;
  //MPI communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  int kk2;

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


  for (l=0;l<bd.NB;l++) {

      if (bd.bnd[l].planex) {
	printf("rank= %d planex %d %d %d\n",rank,bd.bnd[l].ifix[0],bd.bnd[l].ifix[1],bd.bnd[l].ifix[2]);
      }
      else if (bd.bnd[l].planey) {
	printf("rank= %d planey %d %d %d\n",rank,bd.bnd[l].ifix[0],bd.bnd[l].ifix[1],bd.bnd[l].ifix[2]);
      }
      else if (bd.bnd[l].planez) {
	printf("rank= %d planez %d %d %d\n",rank,bd.bnd[l].ifix[0],bd.bnd[l].ifix[1],bd.bnd[l].ifix[2]);
      }
      else if (bd.bnd[l].cylinder) {
	printf("rank= %d cylinder %d %d %d\n",rank,bd.bnd[l].ifix[0],bd.bnd[l].ifix[1],bd.bnd[l].ifix[2]);
      }
      else if (bd.bnd[l].uniform) {
	printf("rank= %d uniform %d %d %d\n",rank,bd.bnd[l].ifix[0],bd.bnd[l].ifix[1],bd.bnd[l].ifix[2]);
      }      
  }


  // loop over all vertices to assign boundary conditions to kfix array
  for ( auto v : verts ) {

    node[v].mp=0.0;
    node[v].un[0]=node[v].un[1]=node[v].un[2]=0.0;
    node[v].kfix[0]=node[v].kfix[1]=node[v].kfix[2]=0.0;

      //loop over all boundary specifications
    for (l=0;l<bd.NB;l++) {

      if (bd.bnd[l].planex) {
	E=mesh.points[v].px[0]-bd.bnd[l].x;
        if (E<0.0) E=-E;
        if (E<1.e-7) {
          for (j=0;j<=2;j++)
	    if (bd.bnd[l].ifix[j]) node[v].kfix[j]=1.0*bd.bnd[l].ifix[j];
	}
      }
      else if (bd.bnd[l].planey) {
	E=mesh.points[v].px[1]-bd.bnd[l].y;
        if (E<0.0) E=-E;
        if (E<1.e-7) {
          for (j=0;j<=2;j++)
	    if (bd.bnd[l].ifix[j]) node[v].kfix[j]=1.0*bd.bnd[l].ifix[j];
	}
      }
      else if (bd.bnd[l].planez) {
	E=mesh.points[v].px[2]-bd.bnd[l].z;
        if (E<0.0) E=-E;
        if (E<1.e-7) {
          for (j=0;j<=2;j++)
	    if (bd.bnd[l].ifix[j]) node[v].kfix[j]=1.0*bd.bnd[l].ifix[j];
	}
      }
      else if (bd.bnd[l].cylinder) {
        rad=sqrt(mesh.points[v].px[0]*mesh.points[v].px[0]+
                 mesh.points[v].px[1]*mesh.points[v].px[1]);

	E=rad-bd.bnd[l].x;
        if (E<0.0) E=-E;
        if (E<1.e-7) {
          for (j=0;j<=2;j++)
	    if (bd.bnd[l].ifix[j]) node[v].kfix[j]=1.0*bd.bnd[l].ifix[j];
	}
      }
      else if (bd.bnd[l].uniform) {
	for (j=0;j<=2;j++)
	  if (bd.bnd[l].ifix[j]) node[v].kfix[j]=1.0*bd.bnd[l].ifix[j];
	
      }
      

    } //end boundary specifications loop

  } // end vertex loop
    

  //store global node ids at each point
  for ( auto v : verts ) {
    node[v].gid=mesh.points[v].gid;

    iz =  node[v].gid % ((gd.LIMX+1)*(gd.LIMY+1));

    ix = iz % (gd.LIMX+1);

    iy = (iz - ix)/(gd.LIMX+1);

    iz = (node[v].gid -iy*(gd.LIMX+1) -ix)/(gd.LIMX+1)*(gd.LIMY+1);

    node[v].ip=ix;
    node[v].jp=iy;
    node[v].kp=iz;
  }
  
  //setup mesh information and sparse region data structures  
  for ( auto c : cs ) {
    
    //list of global node ids for the cell, in the proper order
    j=mesh.elements[c].points[0];
    cell[c].iord1=node[j].gid;
    j=mesh.elements[c].points[1];
    cell[c].iord2=node[j].gid;
    j=mesh.elements[c].points[2];
    cell[c].iord3=node[j].gid;
    j=mesh.elements[c].points[3];
    cell[c].iord4=node[j].gid;    
    j=mesh.elements[c].points[4];
    cell[c].iord5=node[j].gid;
    j=mesh.elements[c].points[5];
    cell[c].iord6=node[j].gid;
    j=mesh.elements[c].points[6];
    cell[c].iord7=node[j].gid;
    j=mesh.elements[c].points[7];
    cell[c].iord8=node[j].gid;

    //cell global id
    cell[c].id=mesh.elements[c].cid;
    
  }  //end cell loop

   if (gd.input_type==2) { //read from a restart dump

     //read time and jedit values
     sprintf(name,"ensight_%05d/data/time.dat",rank);
     out[0]=fopen(name,"r");
     fscanf(out[0],"%d%d%le",&j,&k,&xa);
     gd.jedit=j;
     gd.ncyc=k;
     gd.time=xa;
     fclose(out[0]);

     jedit=gd.jedit-1;

     sprintf(name,"ensight_%05d/data/engold.%05d.vol",rank,jedit);
     out[1]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[1]))!='\n') {
	 k++;
       }
     } 
     sprintf(name,"ensight_%05d/data/engold.%05d.den",rank,jedit);
     out[2]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[2]))!='\n') {
	 k++;
       }
     }
     sprintf(name,"ensight_%05d/data/engold.%05d.ie",rank,jedit);
     out[3]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[3]))!='\n') {
	 k++;
       }
     }  

  
     for (auto c : cs) {//loop over cells

       fscanf(out[1],"%le",&xa);
       cell[c].vol=xa;
       fscanf(out[2],"%le",&xa);
       cell[c].dc=xa;
       fscanf(out[3],"%le",&xa);
       cell[c].ec=xa;

     } //end cell loop

     
     fclose(out[1]);
     fclose(out[2]);
     fclose(out[3]);


     //read node positions
     sprintf(name,"ensight_%05d/data/engold.%05d.geo",rank,jedit);
     out[0]=fopen(name,"r");
     for (j=1;j<=9;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[0]))!='\n') {
	 k++;
       }
     }

     for (auto v : verts) {
       fscanf(out[0],"%le",&xa);
       mesh.points[v].px[0]=xa;
     }

     for (auto v : verts) {
       fscanf(out[0],"%le",&xa);
       mesh.points[v].px[1]=xa;
     }

     for (auto v : verts) {
       fscanf(out[0],"%le",&xa);
       mesh.points[v].px[2]=xa;
     }

     fclose(out[0]);

     // read velocity data
     sprintf(name,"ensight_%05d/data/engold.%05d.vel",rank,jedit);
     out[0]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[0]))!='\n') {
	 k++;
       }
     }


     for (auto v : verts) {
         fscanf(out[0],"%le",&xa);
         node[v].un[0]=xa;
     }

     for (auto v : verts) {
         fscanf(out[0],"%le",&xa);
         node[v].un[1]=xa;
     }

     for (auto v : verts) {
         fscanf(out[0],"%le",&xa);
         node[v].un[2]=xa;
     }

     fclose(out[0]);


 
            
     //set initial conditions in each cell
     for ( auto c : cs ) {
       
       //find vertices and calculate volume and centroid
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
       } //done ordering cell vertices


        //cell centroid

       i=0;
       cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;

       for ( auto v : mesh.vertices_of_cell(c) ) {
	 x[iorder[i]]=mesh.points[v].px[0];
	 y[iorder[i]]=mesh.points[v].px[1];
	 z[iorder[i]]=mesh.points[v].px[2];
	 cell[c].xc[0]+=0.125*x[iorder[i]];
	 cell[c].xc[1]+=0.125*y[iorder[i]];
	 cell[c].xc[2]+=0.125*z[iorder[i]];
	 cn[iorder[i]]=v;
	 i++;
       }

       cell[c].mass=cell[c].dc*cell[c].vol;

       
     } //end cell loop
     
    

  } //done input_type==2
  else if (gd.input_type==3) { //read from a trilink dump that has been processed into Ensight format
    
     if (gd.trilink_spin)
       num=gd.LIMZ*gd.LIMY;
     else
       num=gd.LIMX*gd.LIMY;

     trilink_cell_t *tri_cell;

     tri_cell = new trilink_cell_t[num];

     //create storage for cell-centered momentum
     px=new double[NE];
     py=new double[NE];
     pz=new double[NE];

     sprintf(name,"trilink/data/engold.vol");
     out[1]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[1]))!='\n') {
	 k++;
       }
     } 
     sprintf(name,"trilink/data/engold.den");
     out[2]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[2]))!='\n') {
	 k++;
       }
     }
     sprintf(name,"trilink/data/engold.ie");
     out[3]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[3]))!='\n') {
	 k++;
       }
     }  


     for (j=0;j<num;j++) {
  
       fscanf(out[1],"%le",&xa);
       tri_cell[j].Vc=xa;
       fscanf(out[2],"%le",&xa);
       tri_cell[j].dc=xa;
       fscanf(out[3],"%le",&xa);
       tri_cell[j].ec=xa;

     } //end trilink cell loop

     fclose(out[1]);
     fclose(out[2]);
     fclose(out[3]);


     //read velocity components
     sprintf(name,"trilink/data/engold.vel0");
     out[0]=fopen(name,"r");
     //printf("rank= %d %s\n",rank,name);
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[0]))!='\n') {
	 k++;
       }
     }
     sprintf(name,"trilink/data/engold.vel1");
     out[1]=fopen(name,"r");
     for (j=1;j<=4;j++) { //skip first four lines in file
       k=0;
       while ((ch=(char)fgetc(out[1]))!='\n') {
	 k++;
       }
     }

     for (j=0;j<num;j++) {    
       fscanf(out[0],"%le",&xa);
       tri_cell[j].vel[0]=xa;
       fscanf(out[1],"%le",&xa);
       tri_cell[j].vel[1]=xa;
     }

     fclose(out[0]);
     fclose(out[1]);


     //set initial conditions in each cell
     for ( auto c : cs ) {

       //find vertices and calculate volume and centroid
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
       } //done ordering cell vertices


       auto & xc1 = cell[c].xc; //cell centroid

       i=0;
       xc1[0]=xc1[1]=xc1[2]=0.0;

       for ( auto v : mesh.vertices_of_cell(c) ) {
	 x[iorder[i]]=mesh.points[v].px[0];
	 y[iorder[i]]=mesh.points[v].px[1];
	 z[iorder[i]]=mesh.points[v].px[2];
	 xc1[0]+=0.125*x[iorder[i]];
	 xc1[1]+=0.125*y[iorder[i]];
	 xc1[2]+=0.125*z[iorder[i]];
	 cn[iorder[i]]=v;
	 i++;
       }

       j=c;       
       
       //cell volume
       cell[c].vol= (x[1]*(y[3]*(-z[0] + z[2]) + y[4]*(z[0] - z[5]) + y[0]*(z[2] + z[3] - z[4] - z[5]) + y[6]*(-z[2] + z[5]) + y[5]*(z[0] - z[2] + z[4] - z[6]) + y[2]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[7]*(y[0]*(-z[3] + z[4]) + y[6]*(z[2] + z[3] - z[4] - z[5]) + y[2]*(z[3] - z[6]) + y[3]*(z[0] - z[2] + z[4] - z[6]) + y[5]*(-z[4] + z[6]) + y[4]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[3]*(y[1]*(z[0] - z[2]) + y[7]*(-z[0] + z[2] - z[4] + z[6]) + y[6]*(z[2] - z[7]) + y[2]*(z[0] + z[1] - z[6] - z[7]) + y[4]*(-z[0] + z[7]) + y[0]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[5]*(y[0]*(z[1] - z[4]) + y[7]*(z[4] - z[6]) + y[2]*(-z[1] + z[6]) + y[1]*(-z[0] + z[2] - z[4] + z[6]) + y[4]*(z[0] + z[1] - z[6] - z[7]) + y[6]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[6]*(y[1]*(z[2] - z[5]) + y[7]*(-z[2] - z[3] + z[4] + z[5]) + y[5]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(z[5] - z[7]) + y[3]*(-z[2] + z[7]) + y[2]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[0]*(y[2]*(z[1] - z[3]) + y[7]*(z[3] - z[4]) + y[5]*(-z[1] + z[4]) + y[1]*(-z[2] - z[3] + z[4] + z[5]) + y[3]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[2]*(y[0]*(-z[1] + z[3]) + y[5]*(z[1] - z[6]) + y[1]*(z[0] + z[3] - z[5] - z[6]) + y[7]*(-z[3] + z[6]) + y[6]*(z[1] - z[3] + z[5] - z[7]) + y[3]*(-z[0] - z[1] + z[6] + z[7])) + 
                   x[4]*(y[1]*(-z[0] + z[5]) + y[7]*(z[0] + z[3] - z[5] - z[6]) + y[3]*(z[0] - z[7]) + y[0]*(z[1] - z[3] + z[5] - z[7]) + y[6]*(-z[5] + z[7]) + y[5]*(-z[0] - z[1] + z[6] + z[7])))/12.;

       if (gd.trilink_spin) { //spin case

	 x1=fd.fill[0].x1;
	 x2=fd.fill[0].x2;
	 y1=fd.fill[0].y1;
	 y2=fd.fill[0].y2;
	 z1=fd.fill[0].z1;
	 z2=fd.fill[0].z2;

	 dx=(x2-x1)/gd.LIMX;
	 dy=(y2-y1)/gd.LIMY;
	 dz=(z2-z1)/gd.LIMZ;

	 rad=sqrt(xc1[0]*xc1[0]+xc1[1]*xc1[1]);

	 theta=acos(xc1[0]/rad);

	 ix=(int)ceil((rad-y1)/dy);
	 iz=(int)ceil((xc1[2]-z1)/dz);

	 if (ix==0) ix==1;
	 if (ix>gd.LIMX) ix=gd.LIMX;

	 cell_id=INDEX_FOR_2D(iz,1,gd.LIMZ,ix,1);

	 cell[c].dc=tri_cell[cell_id].dc;
	 cell[c].ec=tri_cell[cell_id].ec;

	 cell[c].mass=cell[c].dc*cell[c].vol;
	   
	 px[c]=cell[c].mass*tri_cell[cell_id].vel[1]*cos(theta);
	 py[c]=cell[c].mass*tri_cell[cell_id].vel[1]*sin(theta);
	 pz[c]=cell[c].mass*tri_cell[cell_id].vel[0];

       }
       else { //extrude case

       
	 //find global i,j,k indices
	 iz = cell[c].id % (gd.LIMX*gd.LIMY);

	 ix = iz % gd.LIMX;

         iy = (iz - ix)/gd.LIMX;

	 iz = (cell[c].id -iy*gd.LIMX -ix)/(gd.LIMX*gd.LIMY);

	 ix++;
	 iy++;
	 iz++;     

	 cell_id=INDEX_FOR_2D(ix,1,gd.LIMX,iy,1);

	 cell[c].dc=tri_cell[cell_id].dc;
	 cell[c].ec=tri_cell[cell_id].ec;

	 cell[c].mass=cell[c].dc*cell[c].vol;	 
       
	 px[c]=cell[c].mass*tri_cell[cell_id].vel[0];
	 py[c]=cell[c].mass*tri_cell[cell_id].vel[1];
	 pz[c]=0.0;

       }
       
     } //end cell loop

     delete[] tri_cell; 

     //find velocities at nodes
     for ( auto v : verts ) {
       node[v].mp=0.0;
       node[v].un[0]=node[v].un[1]=node[v].un[2]=0.0;
     }


     for ( auto v : verts ) { 

       if (vflag[v]==1) //gather corner masses and momenta
	 for ( auto c : mesh.cells_of_point(v) ) {

	   node[v].mp+=0.125*cell[c].mass;
	   node[v].un[0]+=0.125*px[c];
	   node[v].un[1]+=0.125*py[c];
	   node[v].un[2]+=0.125*pz[c];

	 }

     }

     for ( auto v : verts )  
       if (vflag[v]==1) {
	 node[v].un[0]/=node[v].mp;
	 node[v].un[1]/=node[v].mp;
	 node[v].un[2]/=node[v].mp;
       }

     delete[] px;
     delete[] py;
     delete[] pz;


     
  } //done input_type==3  

  else { //use fill commands from the input deck to set initial conditions  
  //loop over all fill instructions in the input deck
  for (l=0;l<fd.NF;l++) {

    if (fd.fill[l].type==4) {

      /* define the spatial extent of the region to fill */
      x1=fd.fill[l].x1;
      x2=fd.fill[l].x2;
      y1=fd.fill[l].y1;
      y2=fd.fill[l].y2;
      z1=fd.fill[l].z1;
      z2=fd.fill[l].z2;
     
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
        cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;

	for ( auto v : mesh.vertices_of_cell(c) ) {
	  x[iorder[i]]=mesh.points[v].px[0];
          y[iorder[i]]=mesh.points[v].px[1];
          z[iorder[i]]=mesh.points[v].px[2];
          cell[c].xc[0]+=0.125*x[iorder[i]];
          cell[c].xc[1]+=0.125*y[iorder[i]];
          cell[c].xc[2]+=0.125*z[iorder[i]];
	  cn[iorder[i]]=v;

	  //    if (rank==0) printf("%d %d %f %f %f\n",c, v,x[iorder[i]],y[iorder[i]],z[iorder[i]]);

	  i++;
	}

	j=c;

        //check if cell is within the fill region
        if (cell[c].xc[0]>=x1 && cell[c].xc[0]<=x2 &&
            cell[c].xc[1]>=y1 && cell[c].xc[1]<=y2 &&
            cell[c].xc[2]>=z1 && cell[c].xc[2]<=z2) {

            cell[c].vol= (x[1]*(y[3]*(-z[0] + z[2]) + y[4]*(z[0] - z[5]) + y[0]*(z[2] + z[3] - z[4] - z[5]) + y[6]*(-z[2] + z[5]) + y[5]*(z[0] - z[2] + z[4] - z[6]) + y[2]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[7]*(y[0]*(-z[3] + z[4]) + y[6]*(z[2] + z[3] - z[4] - z[5]) + y[2]*(z[3] - z[6]) + y[3]*(z[0] - z[2] + z[4] - z[6]) + y[5]*(-z[4] + z[6]) + y[4]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[3]*(y[1]*(z[0] - z[2]) + y[7]*(-z[0] + z[2] - z[4] + z[6]) + y[6]*(z[2] - z[7]) + y[2]*(z[0] + z[1] - z[6] - z[7]) + y[4]*(-z[0] + z[7]) + y[0]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[5]*(y[0]*(z[1] - z[4]) + y[7]*(z[4] - z[6]) + y[2]*(-z[1] + z[6]) + y[1]*(-z[0] + z[2] - z[4] + z[6]) + y[4]*(z[0] + z[1] - z[6] - z[7]) + y[6]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[6]*(y[1]*(z[2] - z[5]) + y[7]*(-z[2] - z[3] + z[4] + z[5]) + y[5]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(z[5] - z[7]) + y[3]*(-z[2] + z[7]) + y[2]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[0]*(y[2]*(z[1] - z[3]) + y[7]*(z[3] - z[4]) + y[5]*(-z[1] + z[4]) + y[1]*(-z[2] - z[3] + z[4] + z[5]) + y[3]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[2]*(y[0]*(-z[1] + z[3]) + y[5]*(z[1] - z[6]) + y[1]*(z[0] + z[3] - z[5] - z[6]) + y[7]*(-z[3] + z[6]) + y[6]*(z[1] - z[3] + z[5] - z[7]) + y[3]*(-z[0] - z[1] + z[6] + z[7])) + 
                   x[4]*(y[1]*(-z[0] + z[5]) + y[7]*(z[0] + z[3] - z[5] - z[6]) + y[3]*(z[0] - z[7]) + y[0]*(z[1] - z[3] + z[5] - z[7]) + y[6]*(-z[5] + z[7]) + y[5]*(-z[0] - z[1] + z[6] + z[7])))/12.;


	    cell[c].mass=fd.fill[l].rho*cell[c].vol;
            cell[c].dc=fd.fill[l].rho;
	    cell[c].ec=fd.fill[l].en;


	}

      } //end cell loop

    } //end if type=4

    else if (fd.fill[l].type==5) { //specify initial velocity in a region

      /* define the spatial extent of the region to fill */
      x1=fd.fill[l].x1;
      x2=fd.fill[l].x2;
      y1=fd.fill[l].y1;
      y2=fd.fill[l].y2;
      z1=fd.fill[l].z1;
      z2=fd.fill[l].z2;

     
      for ( auto v : verts ) {

        
        //check if cell is within the fill region
        if (mesh.points[v].px[0]>=x1 && mesh.points[v].px[0]<=x2 &&
            mesh.points[v].px[1]>=y1 && mesh.points[v].px[1]<=y2 &&
            mesh.points[v].px[2]>=z1 && mesh.points[v].px[2]<=z2) {

	    node[v].un[0]=fd.fill[l].u;
	    node[v].un[1]=fd.fill[l].v;
	    node[v].un[2]=fd.fill[l].w;

	}
      }

    } //end if type=5

    else if (fd.fill[l].type==6) {

      /* noh problem setup */
      /* loop over all vertices on the mesh */
     
      for ( auto v : verts ) {

	if (mesh.points[v].px[0]==0.0 && mesh.points[v].px[1]==0.0 && mesh.points[v].px[2]==0.0) {

	  node[v].un[0]=0.0;
	  node[v].un[1]=0.0;
	  node[v].un[2]=0.0;            

	}
	else {

	  E=sqrt(mesh.points[v].px[0]*mesh.points[v].px[0]+mesh.points[v].px[1]*mesh.points[v].px[1]+mesh.points[v].px[2]*mesh.points[v].px[2]);

	  node[v].un[0]=-1.0*mesh.points[v].px[0]/E;
	  node[v].un[1]=-1.0*mesh.points[v].px[1]/E;
	  node[v].un[2]=-1.0*mesh.points[v].px[2]/E;
	}

      }

    } //end if type =6


    else if (fd.fill[l].type==7) {

      /* verney problem setup */
      /* loop over all vertices on the mesh */

     
      for ( auto v : verts ) {

        
	E=sqrt(mesh.points[v].px[0]*mesh.points[v].px[0]+mesh.points[v].px[1]*mesh.points[v].px[1]+mesh.points[v].px[2]*mesh.points[v].px[2]);

	E=E*E*E;

	node[v].un[0]=fd.fill[l].u*mesh.points[v].px[0]/E;
	node[v].un[1]=fd.fill[l].u*mesh.points[v].px[1]/E;
	node[v].un[2]=fd.fill[l].u*mesh.points[v].px[2]/E;

      }

    } //end type=7

    else if (fd.fill[l].type==8) {

      /* taylor-green vortex problem setup */
      /* loop over all vertices on the mesh */

     
      for ( auto v : verts ) {

	E=sqrt(mesh.points[v].px[0]*mesh.points[v].px[0]+mesh.points[v].px[1]*mesh.points[v].px[1]+mesh.points[v].px[2]*mesh.points[v].px[2]);

	node[v].un[0]=sin(Pi*mesh.points[v].px[0])*cos(Pi*mesh.points[v].px[1]);
	node[v].un[1]=-sin(Pi*mesh.points[v].px[1])*cos(Pi*mesh.points[v].px[0]);
	node[v].un[2]=0.0;

      }

      /* setup region properties */
      /* loop over all cells on rank */

      
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
        
        //find centroid
	i=0;
        cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;

	for ( auto v : mesh.vertices_of_cell(c) ) {
	  x[iorder[i]]=mesh.points[v].px[0];
          y[iorder[i]]=mesh.points[v].px[1];
          z[iorder[i]]=mesh.points[v].px[2];
          cell[c].xc[0]+=0.125*x[iorder[i]];
          cell[c].xc[1]+=0.125*y[iorder[i]];
          cell[c].xc[2]+=0.125*z[iorder[i]];
	  i++;
	}

        cell[c].vol= (x[1]*(y[3]*(-z[0] + z[2]) + y[4]*(z[0] - z[5]) + y[0]*(z[2] + z[3] - z[4] - z[5]) + y[6]*(-z[2] + z[5]) + y[5]*(z[0] - z[2] + z[4] - z[6]) + y[2]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[7]*(y[0]*(-z[3] + z[4]) + y[6]*(z[2] + z[3] - z[4] - z[5]) + y[2]*(z[3] - z[6]) + y[3]*(z[0] - z[2] + z[4] - z[6]) + y[5]*(-z[4] + z[6]) + y[4]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[3]*(y[1]*(z[0] - z[2]) + y[7]*(-z[0] + z[2] - z[4] + z[6]) + y[6]*(z[2] - z[7]) + y[2]*(z[0] + z[1] - z[6] - z[7]) + y[4]*(-z[0] + z[7]) + y[0]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[5]*(y[0]*(z[1] - z[4]) + y[7]*(z[4] - z[6]) + y[2]*(-z[1] + z[6]) + y[1]*(-z[0] + z[2] - z[4] + z[6]) + y[4]*(z[0] + z[1] - z[6] - z[7]) + y[6]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[6]*(y[1]*(z[2] - z[5]) + y[7]*(-z[2] - z[3] + z[4] + z[5]) + y[5]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(z[5] - z[7]) + y[3]*(-z[2] + z[7]) + y[2]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[0]*(y[2]*(z[1] - z[3]) + y[7]*(z[3] - z[4]) + y[5]*(-z[1] + z[4]) + y[1]*(-z[2] - z[3] + z[4] + z[5]) + y[3]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[2]*(y[0]*(-z[1] + z[3]) + y[5]*(z[1] - z[6]) + y[1]*(z[0] + z[3] - z[5] - z[6]) + y[7]*(-z[3] + z[6]) + y[6]*(z[1] - z[3] + z[5] - z[7]) + y[3]*(-z[0] - z[1] + z[6] + z[7])) + 
                   x[4]*(y[1]*(-z[0] + z[5]) + y[7]*(z[0] + z[3] - z[5] - z[6]) + y[3]*(z[0] - z[7]) + y[0]*(z[1] - z[3] + z[5] - z[7]) + y[6]*(-z[5] + z[7]) + y[5]*(-z[0] - z[1] + z[6] + z[7])))/12.;


	cell[c].mass=cell[c].vol;//density is 1g per cc
	cell[c].dc=1.0;
	cell[c].ec=(10.+0.25*(cos(2.*Pi*cell[c].xc[0])+cos(2.*Pi*cell[c].xc[1])))/(2./3.);


      } //end cell loop

    } //end type=8

    else if (fd.fill[l].type==9) {

      /* kidder problem setup */
      /* loop over all vertices on the mesh */
      for ( auto v : verts ) {

	E=sqrt(mesh.points[v].px[0]*mesh.points[v].px[0]+mesh.points[v].px[1]*mesh.points[v].px[1]+mesh.points[v].px[2]*mesh.points[v].px[2]);

	node[v].un[0]=-0.5*mesh.points[v].px[0];
	node[v].un[1]=-0.5*mesh.points[v].px[1];
	node[v].un[2]=-0.5*mesh.points[v].px[2];	

      }      
     

      /* setup region properties */
     /* loop over all cells on rank */
 
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

        //find centroid
	i=0;
        cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;

	for ( auto v : mesh.vertices_of_cell(c) ) {
	  x[iorder[i]]=mesh.points[v].px[0];
          y[iorder[i]]=mesh.points[v].px[1];
          z[iorder[i]]=mesh.points[v].px[2];
          cell[c].xc[0]+=0.125*x[iorder[i]];
          cell[c].xc[1]+=0.125*y[iorder[i]];
          cell[c].xc[2]+=0.125*z[iorder[i]];
	  i++;
	}
	  
	cell[c].vol= (x[1]*(y[3]*(-z[0] + z[2]) + y[4]*(z[0] - z[5]) + y[0]*(z[2] + z[3] - z[4] - z[5]) + y[6]*(-z[2] + z[5]) + y[5]*(z[0] - z[2] + z[4] - z[6]) + y[2]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[7]*(y[0]*(-z[3] + z[4]) + y[6]*(z[2] + z[3] - z[4] - z[5]) + y[2]*(z[3] - z[6]) + y[3]*(z[0] - z[2] + z[4] - z[6]) + y[5]*(-z[4] + z[6]) + y[4]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[3]*(y[1]*(z[0] - z[2]) + y[7]*(-z[0] + z[2] - z[4] + z[6]) + y[6]*(z[2] - z[7]) + y[2]*(z[0] + z[1] - z[6] - z[7]) + y[4]*(-z[0] + z[7]) + y[0]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[5]*(y[0]*(z[1] - z[4]) + y[7]*(z[4] - z[6]) + y[2]*(-z[1] + z[6]) + y[1]*(-z[0] + z[2] - z[4] + z[6]) + y[4]*(z[0] + z[1] - z[6] - z[7]) + y[6]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[6]*(y[1]*(z[2] - z[5]) + y[7]*(-z[2] - z[3] + z[4] + z[5]) + y[5]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(z[5] - z[7]) + y[3]*(-z[2] + z[7]) + y[2]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[0]*(y[2]*(z[1] - z[3]) + y[7]*(z[3] - z[4]) + y[5]*(-z[1] + z[4]) + y[1]*(-z[2] - z[3] + z[4] + z[5]) + y[3]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[2]*(y[0]*(-z[1] + z[3]) + y[5]*(z[1] - z[6]) + y[1]*(z[0] + z[3] - z[5] - z[6]) + y[7]*(-z[3] + z[6]) + y[6]*(z[1] - z[3] + z[5] - z[7]) + y[3]*(-z[0] - z[1] + z[6] + z[7])) +
		      x[4]*(y[1]*(-z[0] + z[5]) + y[7]*(z[0] + z[3] - z[5] - z[6]) + y[3]*(z[0] - z[7]) + y[0]*(z[1] - z[3] + z[5] - z[7]) + y[6]*(-z[5] + z[7]) + y[5]*(-z[0] - z[1] + z[6] + z[7])))/12.;

	cell[c].dc=exp(-0.5*(cell[c].xc[0]*cell[c].xc[0]+cell[c].xc[1]*cell[c].xc[1]+cell[c].xc[2]*cell[c].xc[2]))/sqrt(2.);
	cell[c].mass=cell[c].dc*cell[c].vol;
	cell[c].ec=3./8.;
	    

      } //end cell loop

    } //end type=9
 
    else if (fd.fill[l].type==10) { //setup of cylindrical region

      /* define the spatial extent of the region to fill */
      x1=fd.fill[l].theta1*Pi/180.0;
      x2=fd.fill[l].theta2*Pi/180.0;
      y1=fd.fill[l].rad1;
      y2=fd.fill[l].rad2;
      z1=fd.fill[l].z1;
      z2=fd.fill[l].z2;

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
        cell[c].xc[0]=cell[c].xc[1]=cell[c].xc[2]=0.0;

	for ( auto v : mesh.vertices_of_cell(c) ) {
	  x[iorder[i]]=mesh.points[v].px[0];
          y[iorder[i]]=mesh.points[v].px[1];
          z[iorder[i]]=mesh.points[v].px[2];
          cell[c].xc[0]+=0.125*x[iorder[i]];
          cell[c].xc[1]+=0.125*y[iorder[i]];
          cell[c].xc[2]+=0.125*z[iorder[i]];
	  cn[iorder[i]]=v;

	  i++;
	}

	j=c;
        rad=sqrt(cell[c].xc[0]*cell[c].xc[0]+cell[c].xc[1]*cell[c].xc[1]);
        theta1=acos(cell[c].xc[0]/rad);

        //check if cell is within the fill region
        if (theta1<=x1 && theta1>=x2 &&
            rad>=y1 && rad<=y2 &&
            cell[c].xc[2]>=z1 && cell[c].xc[2]<=z2) {

            cell[c].vol= (x[1]*(y[3]*(-z[0] + z[2]) + y[4]*(z[0] - z[5]) + y[0]*(z[2] + z[3] - z[4] - z[5]) + y[6]*(-z[2] + z[5]) + y[5]*(z[0] - z[2] + z[4] - z[6]) + y[2]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[7]*(y[0]*(-z[3] + z[4]) + y[6]*(z[2] + z[3] - z[4] - z[5]) + y[2]*(z[3] - z[6]) + y[3]*(z[0] - z[2] + z[4] - z[6]) + y[5]*(-z[4] + z[6]) + y[4]*(-z[0] - z[3] + z[5] + z[6])) + 
                   x[3]*(y[1]*(z[0] - z[2]) + y[7]*(-z[0] + z[2] - z[4] + z[6]) + y[6]*(z[2] - z[7]) + y[2]*(z[0] + z[1] - z[6] - z[7]) + y[4]*(-z[0] + z[7]) + y[0]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[5]*(y[0]*(z[1] - z[4]) + y[7]*(z[4] - z[6]) + y[2]*(-z[1] + z[6]) + y[1]*(-z[0] + z[2] - z[4] + z[6]) + y[4]*(z[0] + z[1] - z[6] - z[7]) + y[6]*(-z[1] - z[2] + z[4] + z[7])) + 
                   x[6]*(y[1]*(z[2] - z[5]) + y[7]*(-z[2] - z[3] + z[4] + z[5]) + y[5]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(z[5] - z[7]) + y[3]*(-z[2] + z[7]) + y[2]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[0]*(y[2]*(z[1] - z[3]) + y[7]*(z[3] - z[4]) + y[5]*(-z[1] + z[4]) + y[1]*(-z[2] - z[3] + z[4] + z[5]) + y[3]*(z[1] + z[2] - z[4] - z[7]) + y[4]*(-z[1] + z[3] - z[5] + z[7])) + 
                   x[2]*(y[0]*(-z[1] + z[3]) + y[5]*(z[1] - z[6]) + y[1]*(z[0] + z[3] - z[5] - z[6]) + y[7]*(-z[3] + z[6]) + y[6]*(z[1] - z[3] + z[5] - z[7]) + y[3]*(-z[0] - z[1] + z[6] + z[7])) + 
                   x[4]*(y[1]*(-z[0] + z[5]) + y[7]*(z[0] + z[3] - z[5] - z[6]) + y[3]*(z[0] - z[7]) + y[0]*(z[1] - z[3] + z[5] - z[7]) + y[6]*(-z[5] + z[7]) + y[5]*(-z[0] - z[1] + z[6] + z[7])))/12.;

	    cell[c].mass=fd.fill[l].rho*cell[c].vol;
            cell[c].dc=fd.fill[l].rho;
	    cell[c].ec=fd.fill[l].en;

	} // end within region test

      } //end cell loop

    } //end if type=10
 

  } /* end fill type loop */

  } //end reading setup instructions from input deck



  for ( auto c : cs ) {

    cell[c].flg=0.0;
    for (i=0;i<3;i++) {
      cell[c].stress1[i]=0.0;
      cell[c].stress2[i]=0.0;
      cell[c].stress3[i]=0.0;
    }
      
  }


  for ( auto c : cs_owned ) {

    cell[c].flg=1.0;

  }


/* calculate the vertex masses */
/* loop over all vertices on rank */
  for ( auto v : verts ) 
    node[v].mp=0.0;

  for ( auto v : verts ) {  

    if (vflag[v]==1) //gather corner masses
      for ( auto c : mesh.cells_of_point(v) ) {

	node[v].mp+=0.125*cell[c].mass;

      }

  }

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

	   //nodal masses and velocities
	   for (auto cr : mesh.corners_of_cell(c)) {
	     
	     auto v = mesh.corners[cr].point;
	      
	     send[ndat*i+index+offset]=node[v].mp;
	     index++;

	     send[ndat*i+index+offset]=node[v].un[0];
	     index++;

	     send[ndat*i+index+offset]=node[v].un[1];
	     index++;

	     send[ndat*i+index+offset]=node[v].un[2];
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

	   //nodal masses and velocities
	   for (int cr=0; cr<8; cr++) {
	     
	     auto v = verts[cn[cr]];
	       
	     //only update vertex position for ghost nodes
	     if (vflag[v]==0) node[v].mp=recv[ndat*i+index+offset];
	     index++;
	     if (vflag[v]==0) node[v].un[0]=recv[ndat*i+index+offset];
	     index++;
	     if (vflag[v]==0) node[v].un[1]=recv[ndat*i+index+offset];
	     index++;
	     if (vflag[v]==0) node[v].un[2]=recv[ndat*i+index+offset];
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
