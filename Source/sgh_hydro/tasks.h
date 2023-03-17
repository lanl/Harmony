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
#define INDEX_FOR_3D(x,xl,xh,y,yl,yh,z,zl) (((xh)-(xl)+1)*((y)-(yl))+(x)-(xl) + ((xh)-(xl)+1)*((yh)-(yl)+1)*((z)-(zl)))
#define INDEX_FOR_2D(x,xl,xh,y,yl) (((xh)-(xl)+1)*((y)-(yl))+(x)-(xl))
#define tau(a,b) tau[INDEX_FOR_2D(a,1,3,b,1)]
#define D(a,b) D[INDEX_FOR_2D(a,1,3,b,1)]
#define W(a,b) W[INDEX_FOR_2D(a,1,3,b,1)]
#define areab(a,b) areab[INDEX_FOR_2D(a,1,3,b,1)]
#define ac(a,b) ac[INDEX_FOR_2D(a,1,3,b,1)]

#define MAX(a,b) ((b)>(a) ? (b) : (a))
#define MIN(a,b) ((b)>(a) ? (a) : (b))
#define SIGN(a,b) ((b)==0. ? 0. : (b)/fabs(b))
#define COPYSIGN(a,b) (fabs(a)*(b)/fabs(b))

#include "mesh_type.h"
#include "types.h"


typedef struct _contour {
  char name[50];
  int points;
} contour_t;

typedef struct _contour_data {
  contour_t contour[20];
  int NC;
} con_data_t;

typedef struct _fill {
  int type;
  int j1,j2;
  int i1,i2;
  int k1,k2;
  double x1,x2;
  double y1,y2;
  double z1,z2;
  int con1,con2;
  double rad1,rad2;
  double theta1,theta2;
  double phi1,phi2;
  int reg1,reg2;
  double u,v,w;
  double en,rho;
  int num_pts;
  std::vector<double> xp,yp;
} fill_t;

typedef struct _fill_data {
  fill_t fill[40];
  int NF;
} fill_data_t;


typedef struct _bound {
  int ifix[3];
  int planex,planey,planez;
  int cylinder;
  int uniform;
  double x,y,z;
} bound_t;

typedef struct _bound_data {

  int NB;
  bound_t bnd[40];

} bound_data_t;

typedef struct _rlx_data { //data for mesh relaxation specifications
  
  int type;
  int nlag;         //number of logical j lines that are specified to be Lagrangian
  int jlines[15];   //list of logical j lines that are specified to be Lagrangian

} rlx_data_t;


typedef struct _global_data {
  int iorder;
  int ihour;
  int edit_type;
  int iremap;
  int isource;
  double TFINAL;
  double DTMIN;
  double DT;
  int igeo;
  double freeze_time;
  double relax_time;
  double etime;
  double CA;
  double C1;
  double C2;
  int LIMX;
  int LIMY;
  int LIMZ;
  int NPEX;
  int NPEY;
  int NPEZ;    
  int max_steps;
  int input_type;
  int kefix;
  int ihgforce;
  
  double time;
  int ifstop;
  int jedit;
  int ncyc;
  int trilink_spin;    
} global_data_t;

typedef struct _mymesh { /* local mesh consisting of 27 cells used to find cell-centered gradients*/

  int indx[27];
  int cpy[27];
  double xc[27],yc[27],zc[27];

} mymesh_t;  


double gas1(mat_t & mat,double d, double t)
{
  double en;

  if (t<0.0)
    en=0.0;
  else
    en=t;

  /* return pressure */
  return (mat.g-1.)*en*d;
}

double gas2(mat_t & mat,double d, double t)
{
  double en,F;

  if (t<0.0)
    en=0.0;
  else
    en=t;

  /* return sound speed */ 
  F=mat.g*(mat.g-1.)*en;

  if (F<mat.csmin) F=mat.csmin;

  return sqrt(F);
}


#include "read_connect.h"
#include "parser.h"
#include "setup.h"
#include "editor.h"
#include "metrics.h"
#include "properties.h"
#include "time.h"
#include "force1.h"
#include "force2.h"
#include "update1.h"
#include "update2.h"

#include "ale_init.h"
#include "ale_init2.h"
#include "complete_stencil.h"
#include "slopes1.h"
#include "slopes2.h"
#include "fluxes1.h"
#include "fluxes2.h"
#include "ale_update.h"
#include "ale_update2.h"
#include "kefix.h"
#include "ale_displace.h"
#include "ale_displace2.h"
#include "rlx_smooth.h"
#include "rlx_smooth2.h"
