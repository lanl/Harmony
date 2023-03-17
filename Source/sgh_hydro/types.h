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
//define class for eos variables and regions of the mesh
class mat_t {

 public:
  
  double cv;
  double g;
  double csmin;
  double pmin;
  
  mat_t();
};

mat_t::mat_t()
{

  cv=1.0;
  g=1.4;
  csmin=0.1;
  pmin=0.0;
  
}



//define class for cell-centered variables 
class cell_t {
  
public :


  // bulk state
  double dc;
  double pc;
  double ec;
  double sound;
  double mass;
  double vol;
  int id;
  double h;
  double phi;
  double flg;
  double de;

  int iord1;
  int iord2;
  int iord3;
  int iord4;
  int iord5;
  int iord6;
  int iord7;
  int iord8;  
  
  double xc[3];
  
  double stress1[3];
  double stress2[3];
  double stress3[3];
  
  //corner forces
  double f1[3];
  double f2[3];
  double f3[3];
  double f4[3];  
  double f5[3];
  double f6[3];
  double f7[3];
  double f8[3];
  
  //b matrix
  double b1[3];
  double b2[3];
  double b3[3];
  double b4[3];  
  double b5[3];
  double b6[3];
  double b7[3];
  double b8[3];  
  
  //gradients  
  double dudx1[3];
  double dudx2[3];
  double dudx3[3];

  
  // ALE variables
  double xcp[3];
  double enc;

  double MFx[3];
  double EFx[3];
  double TFx[3];

  double sx_d[3];
  double sx_ie[3];
  double sx_en[3];

  cell_t();

};

cell_t::cell_t()
{

  dc=0.0;
  pc=0.0;
  ec=0.0;
  sound=0.0;
  mass=0.0;
  vol=0.0;
  id=0;
  h=0.0;
  phi=0.0;
  flg=0.0;
  de=0.0;
  enc=0.0;
  
  iord1=0;
  iord2=0;
  iord3=0;
  iord4=0;
  iord5=0;
  iord6=0;
  iord7=0;
  iord8=0;
  
  for (int i=0;i<3;i++) {
    xc[i]=0.0;
    xcp[i]=0.0;
    
    stress1[i]=0.0;
    stress2[i]=0.0;
    stress3[i]=0.0;
  
    //corner forces
    f1[i]=0.0;
    f2[i]=0.0;
    f3[i]=0.0;
    f4[i]=0.0;
    f5[i]=0.0;
    f6[i]=0.0;
    f7[i]=0.0;
    f8[i]=0.0;
  
    //b matrix
    b1[i]=0.0;
    b2[i]=0.0;
    b3[i]=0.0;
    b4[i]=0.0;
    b5[i]=0.0;
    b6[i]=0.0;
    b7[i]=0.0;
    b8[i]=0.0;
  
    //gradients  
    dudx1[i]=0.0;
    dudx2[i]=0.0;
    dudx3[i]=0.0;

    //ALE variables
    MFx[i]=0.0;
    EFx[i]=0.0;
    TFx[i]=0.0;

    sx_d[i]=0.0;
    sx_ie[i]=0.0;
    sx_en[i]=0.0;     

  }

  

} //end cell constructor




class vertex_t {
  
public :

  double un[3];
  double du[3];
  double dx[3];
  double x0[3];
  double fx[3];
  double mp;
  double kfix[3];
  int gid;
  int ip;
  int jp;
  int kp;
  
  //ALE variables
  double xp[3];
  double MUFx[3];
  double MVFx[3];  
  double MWFx[3];
  double sx_u[3];
  double sx_v[3];  
  double sx_w[3];
  
  vertex_t();

};

vertex_t::vertex_t()
{

  ip=jp=kp=0;
  for (int i=0;i<3;i++) {
    un[i]=0.0;
    du[i]=0.0;
    dx[i]=0.0;
    x0[i]=0.0;
    fx[i]=0.0;
    kfix[i]=0.0;

    xp[i]=0.0;
    MUFx[i]=0.0;
    MVFx[i]=0.0;  
    MWFx[i]=0.0;
    sx_u[i]=0.0;
    sx_v[i]=0.0;  
    sx_w[i]=0.0;    
  }

  mp=0.0;
  gid=0;

}

class mesh_connect_t {
  
public :
  
  std::vector<int> recv_local;
  std::vector<int> recv_rank;  
  std::vector<int> recv_cell;
  std::vector<int> recv_num;  

  std::vector<int> send_local;
  std::vector<int> send_rank;  
  std::vector<int> send_cell;
  std::vector<int> send_num;  

  mesh_connect_t();

};

mesh_connect_t::mesh_connect_t()
{
  //setup initial storage capacity for rank connectivities
  recv_local.reserve(10);
  recv_rank.reserve(10);  
  recv_cell.reserve(10);
  recv_num.reserve(10);
  
  send_local.reserve(10);
  send_rank.reserve(10);  
  send_cell.reserve(10);  
  send_num.reserve(10);  
  

}
