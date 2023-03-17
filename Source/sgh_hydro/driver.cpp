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
#include "tasks.h"

int main(int argc, char** argv) 
{

  // make region data structure
  mat_t mat;

  // make contour data structure
  con_data_t cd;

  // make fill instruction data structure
  fill_data_t fd;

  // make global data structure
  global_data_t gd;

  // make boundary instruction data structure
  bound_data_t bd;

  //ALE mesh relaxation data structure
  rlx_data_t rlx;

  //rank connectivity information
  mesh_connect_t mcon;

  //===========================================================================
  // START MPI
  //===========================================================================  
  int numtasks, rank;
  MPI_Status Stat;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //===========================================================================
  // CREATE MESH
  //===========================================================================  
  char file_name[100];
  sprintf(file_name,"mesh/mesh_%05d.dat",rank);
  
  //initialize the mesh using data from a mesh definition file specific to this rank
  mesh_t mesh(rank);

  MPI_Barrier(MPI_COMM_WORLD);
  
  //default global values
  gd.iorder=0;
  gd.ihour=0;
  gd.edit_type=0;
  gd.iremap=0;
  gd.isource=0;
  gd.TFINAL=5.01;
  gd.DTMIN=1.e-3;
  gd.DT=1.e-3;
  gd.igeo=0;
  gd.freeze_time=100.0;
  gd.relax_time=100.0;
  gd.etime=0.25;
  gd.CA=0.0;
  gd.C1=1.3;
  gd.C2=0.01;
  gd.LIMX=10;
  gd.LIMY=10;
  gd.LIMZ=10;
  gd.max_steps=200000;
  gd.input_type=0;
  gd.kefix=0;
  gd.ihgforce=0;
  
  gd.time=0.0;
  gd.ifstop=0;
  gd.jedit=0;
  gd.ncyc=0;
  gd.trilink_spin=0;  
      
  //default mesh relaxation values
  rlx.type=0; 
  rlx.nlag=0;

  for (int i=0;i<40;i++) {
    fd.fill[i].xp.reserve(1);
    fd.fill[i].yp.reserve(1);
  }
  
  //===========================================================================
  // READ INPUT DECK
  //===========================================================================

  parser(rank,gd,mat,cd,fd,bd,rlx);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("done with parser rank=%d\n",rank);
  
  //===========================================================================
  // Initial conditions
  //===========================================================================
  // the solution time starts at zero
  double soln_time(0);  
  int time_cnt(0); 

  double edtimes[250],t1;
  int num_steps;
  int ifstop=0;

  int material_id,num;

  t1=gd.etime;
  num_steps=0;

  int RKstep,RKorder;
  double RKcoeff;
  
  //create memory storage for cells, nodes, and corners
  std::vector<cell_t> cell;
  std::vector<vertex_t> node;

  num=mesh.elements.size();
  cell.reserve(num);
  cell.resize(num);

  printf("cell number = %d rank=%d\n",num,rank);

  num=mesh.points.size(); 
  node.reserve(num); 
  node.resize(num);

  printf("node number = %d rank=%d\n",num,rank);
  
  printf("ready for setup rank=%d\n",rank);

  //read rank connectivity file for this rank
  read_connect(mcon); 
  
  //setup all quantities on cells and nodes
  setup(mesh,mat,gd,fd,bd,cell,node,mcon);

  MPI_Barrier(MPI_COMM_WORLD);
  time_cnt=num_steps=gd.ncyc;
  
  //find pressure and sound speed for each cell
  properties(mesh,mat,cell);
  
  //time zero edit
  if (gd.jedit==0) {
    if (rank==0) printf("Wrote edit %d at cycle %d and time %e\n",gd.jedit,time_cnt,gd.time);
  
    editor(mesh,mat,cell,node,gd.jedit,time_cnt,gd.time,edtimes);

    gd.jedit++;

    t1=gd.etime;
  }
  else { //start from restart dump

    t1=gd.jedit*gd.etime;

  }

  
  // main physics loop
  for ( ; 
    (num_steps < gd.max_steps && gd.time < gd.TFINAL); 
    ++num_steps 
  ) {


    //calculate vol and b matrix for each cell     
    metrics(mesh,cell,node);

    properties(mesh,mat,cell);
        
    // find time step
    auto local_time_step=time(mesh,cell,node,gd.DT,gd.DTMIN);

    double global_time_step=local_time_step;

    MPI_Allreduce (&local_time_step,&global_time_step,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);       
    gd.DT=global_time_step;
      
    if (rank==0) printf("step =  %d, time = %e, dt = %e\n",time_cnt,gd.time,gd.DT);

    if (gd.time > t1)	{ //do graphics dump

      if (rank==0) printf("Wrote edit %d at cycle %d and time %e\n",gd.jedit,time_cnt,gd.time);
	 
      editor(mesh,mat,cell,node,gd.jedit,time_cnt,gd.time,edtimes);

      gd.jedit++;

      t1+=gd.etime;
    }

    //predictor step
    force1(mesh,mat,cell,node,mcon,gd.DT,gd.C1,gd.C2,gd.CA);

    update1(mesh,mat,gd.time,gd.DT,time_cnt,cell,node,mcon);

    //corrector step
    metrics(mesh,cell,node);

    properties(mesh,mat,cell);
    
    force2(mesh,mat,cell,node,mcon,gd.DT,gd.C1,gd.C2,gd.CA,gd.ihgforce);

    update2(mesh,mat,gd.time,gd.DT,time_cnt,cell,node,mcon);
    

    //remap step
     if (gd.iremap) {

    
       ale_init(mesh,node);
		
       if (gd.iremap==2) { //Winslow mesh relaxation

	 ale_displace2(mesh,node);

	 for (int l=1;l<=5;l++) { // 5 iterations

	   if (rlx.type==1) { //ALE mesh relaxation on cylindrical mesh

	     rlx_smooth2(mesh,gd,rlx,cell,node);
	   }
	   else { //default ALE mesh relaxation on box mesh for Noh type problem

	     rlx_smooth(mesh,gd,cell,node);
	   }

	 } //end mesh smoothing iterations

       } //end Winlsow mesh relaxation
 
     
       ale_init2(mesh,mat,cell,node);

      //calculation gradients for higher order remap
       if (gd.iorder) {
    
	 slopes1(mesh,mat,gd,cell,node,mcon);

	 slopes2(mesh,cell,node,mcon);

       }
       
       //calculate fluxes at each face
       fluxes1(mesh,mat,cell,node,mcon);

       fluxes2(mesh,cell,node,mcon);
       
       //update vertex positions
       ale_displace(mesh,node);
       
       //update material mass momentum and energy in the cells
       ale_update(mesh,mat,cell,node,mcon);

       ale_update2(mesh,cell,node,mcon);

       //use kefix to find new internal energy
       if (gd.kefix) kefix(mesh,mat,cell,node,mcon);
       
    } //end remap step
    
   
   
    // update time
    gd.time = gd.time + gd.DT;
    time_cnt = time_cnt+1;

  } //end main physics loop
  
  MPI_Finalize();

  return 0;  

}

