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
void rlx_smooth2( 
  mesh_t  &mesh,
  global_data_t &gd,
  rlx_data_t &rlx,
  std::vector<cell_t> &cell,
  std::vector<vertex_t> &node   
) {

  /* smooths vertices according to Winslow method for cylindrical problems */

  double factor,phi[4],fmin[4],fmax[4],xc[8],yc[8],zc[8];

  int m,n,l,i,j,k;

  int iorder[8],cn[8],myvert[8];

  int indx[27];


  double a1,a2,a3,b1,b2,b3;

  double rxx,ryy,rzz,rxy,rxz,ryz;

  double ux,uy,uz,uxx,uyy,uzz,uxy,uxz,uyz;

  double vx,vy,vz;

  double wx,wy,wz;

  // get the context
  int rank=0;

  //cells
  auto cs = mesh.cells();
  auto verts = mesh.vertices();
  auto v1 = verts[0];

  /* loop over all verticies */
  for ( auto v : mesh.vertices() ) {

      for (j=0;j<3;j++) {
	node[v].x0[j]=node[v].xp[j];
      }

  }

  /* loop over all verticies */
  for ( auto v : mesh.vertices() ) {

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
    // now all entries in the indx array that are not -1 have been set
 
    // begin smoothing the vertices, starting with the interior vertices
    if (node[v].ip!=0 && node[v].ip!=gd.LIMX && node[v].kp!=0 && node[v].kp!=gd.LIMZ && node[v].jp!=0 && node[v].jp!=gd.LIMY) {  

      ux=(node[indx[14]].xp[0]-node[indx[12]].xp[0]);
      uy=(node[indx[16]].xp[0]-node[indx[10]].xp[0]);
      uz=(node[indx[22]].xp[0]-node[indx[4]].xp[0]);

      vx=(node[indx[14]].xp[1]-node[indx[12]].xp[1]);
      vy=(node[indx[16]].xp[1]-node[indx[10]].xp[1]);
      vz=(node[indx[22]].xp[1]-node[indx[4]].xp[1]);

      wx=(node[indx[14]].xp[2]-node[indx[12]].xp[2]);
      wy=(node[indx[16]].xp[2]-node[indx[10]].xp[2]);
      wz=(node[indx[22]].xp[2]-node[indx[4]].xp[2]);

      rxx=ux*ux+vx*vx+wx*wx;
      ryy=uy*uy+vy*vy+wy*wy;
      rzz=uz*uz+vz*vz+wz*wz;

      rxy=ux*uy+vx*vy+wx*wy;
      rxz=ux*uz+vx*vz+wx*wz;
      ryz=uz*uy+vz*vy+wz*wy;

      a1=ryy*rzz-ryz*ryz;
      a2=rzz*rxx-rxz*rxz;
      a3=rxx*ryy-rxy*rxy;

      b1=0.5*(ryz*rxz-rxy*rzz);
      b2=0.5*(rxz*rxy-ryz*rxx);
      b3=0.5*(rxy*ryz-rxz*ryy);

      uxx=(node[indx[14]].xp[0]+node[indx[12]].xp[0]);
      uyy=(node[indx[16]].xp[0]+node[indx[10]].xp[0]);
      uzz=(node[indx[22]].xp[0]+node[indx[4]].xp[0]);

      uxy=(node[indx[17]].xp[0]-node[indx[15]].xp[0]+node[indx[9]].xp[0]-node[indx[11]].xp[0]);
      uyz=(node[indx[26]].xp[0]-node[indx[19]].xp[0]+node[indx[1]].xp[0]-node[indx[7]].xp[0]);
      uxz=(node[indx[23]].xp[0]-node[indx[21]].xp[0]+node[indx[3]].xp[0]-node[indx[5]].xp[0]);

 
      node[v].x0[0]=(a1*uxx+a2*uyy+a3*uzz+b1*uxy+b2*uyz+b3*uxz)
	/(2.0*(a1+a2+a3));

      uxx=(node[indx[14]].xp[1]+node[indx[12]].xp[1]);
      uyy=(node[indx[16]].xp[1]+node[indx[10]].xp[1]);
      uzz=(node[indx[22]].xp[1]+node[indx[4]].xp[1]);

      uxy=(node[indx[17]].xp[1]-node[indx[15]].xp[1]+node[indx[9]].xp[1]-node[indx[11]].xp[1]);
      uyz=(node[indx[26]].xp[1]-node[indx[19]].xp[1]+node[indx[1]].xp[1]-node[indx[7]].xp[1]);
      uxz=(node[indx[23]].xp[1]-node[indx[21]].xp[1]+node[indx[3]].xp[1]-node[indx[5]].xp[1]);

      node[v].x0[1]=(a1*uxx+a2*uyy+a3*uzz+b1*uxy+b2*uyz+b3*uxz)
	/(2.0*(a1+a2+a3));

      uxx=(node[indx[14]].xp[2]+node[indx[12]].xp[2]);
      uyy=(node[indx[16]].xp[2]+node[indx[10]].xp[2]);
      uzz=(node[indx[22]].xp[2]+node[indx[4]].xp[2]);

      uxy=(node[indx[17]].xp[2]-node[indx[15]].xp[2]+node[indx[9]].xp[2]-node[indx[11]].xp[2]);
      uyz=(node[indx[26]].xp[2]-node[indx[19]].xp[2]+node[indx[1]].xp[2]-node[indx[7]].xp[2]);
      uxz=(node[indx[23]].xp[2]-node[indx[21]].xp[2]+node[indx[3]].xp[2]-node[indx[5]].xp[2]);

      node[v].x0[2]=(a1*uxx+a2*uyy+a3*uzz+b1*uxy+b2*uyz+b3*uxz)
          	/(2.0*(a1+a2+a3));

    }
    else if (node[v].ip==0 && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp!=0 && node[v].kp!=gd.LIMZ) { // -x plane
      node[v].x0[0]=0.125*(node[indx[1]].xp[0]+node[indx[10]].xp[0]+node[indx[19]].xp[0]+node[indx[22]].xp[0]+
		      node[indx[25]].xp[0]+node[indx[16]].xp[0]+node[indx[7]].xp[0]+node[indx[4]].xp[0]);
      node[v].x0[1]=0.125*(node[indx[1]].xp[1]+node[indx[10]].xp[1]+node[indx[19]].xp[1]+node[indx[22]].xp[1]+
		      node[indx[25]].xp[1]+node[indx[16]].xp[1]+node[indx[7]].xp[1]+node[indx[4]].xp[1]);
      node[v].x0[2]=0.125*(node[indx[1]].xp[2]+node[indx[10]].xp[2]+node[indx[19]].xp[2]+node[indx[22]].xp[2]+
		      node[indx[25]].xp[2]+node[indx[16]].xp[2]+node[indx[7]].xp[2]+node[indx[4]].xp[2]);
    }

    else if (node[v].ip==gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp!=0 && node[v].kp!=gd.LIMZ) { // +x plane
      node[v].x0[0]=0.125*(node[indx[1]].xp[0]+node[indx[10]].xp[0]+node[indx[19]].xp[0]+node[indx[22]].xp[0]+
		      node[indx[25]].xp[0]+node[indx[16]].xp[0]+node[indx[7]].xp[0]+node[indx[4]].xp[0]);
      node[v].x0[1]=0.125*(node[indx[1]].xp[1]+node[indx[10]].xp[1]+node[indx[19]].xp[1]+node[indx[22]].xp[1]+
		      node[indx[25]].xp[1]+node[indx[16]].xp[1]+node[indx[7]].xp[1]+node[indx[4]].xp[1]);
      node[v].x0[2]=0.125*(node[indx[1]].xp[2]+node[indx[10]].xp[2]+node[indx[19]].xp[2]+node[indx[22]].xp[2]+
		      node[indx[25]].xp[2]+node[indx[16]].xp[2]+node[indx[7]].xp[2]+node[indx[4]].xp[2]);
    }

    else if (node[v].ip!=0 && node[v].ip!=gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp==0 && node[v].kp!=gd.LIMZ) { // -z plane
      node[v].x0[0]=0.125*(node[indx[11]].xp[0]+node[indx[14]].xp[0]+node[indx[17]].xp[0]+node[indx[16]].xp[0]+
		      node[indx[15]].xp[0]+node[indx[12]].xp[0]+node[indx[9]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.125*(node[indx[11]].xp[1]+node[indx[14]].xp[1]+node[indx[17]].xp[1]+node[indx[16]].xp[1]+
		      node[indx[15]].xp[1]+node[indx[12]].xp[1]+node[indx[9]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.125*(node[indx[11]].xp[2]+node[indx[14]].xp[2]+node[indx[17]].xp[2]+node[indx[16]].xp[2]+
		      node[indx[15]].xp[2]+node[indx[12]].xp[2]+node[indx[9]].xp[2]+node[indx[10]].xp[2]);
    }

    else if (node[v].ip!=0 && node[v].ip!=gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp!=0 && node[v].kp==gd.LIMZ) { // +z plane
      node[v].x0[0]=0.125*(node[indx[11]].xp[0]+node[indx[14]].xp[0]+node[indx[17]].xp[0]+node[indx[16]].xp[0]+
		      node[indx[15]].xp[0]+node[indx[12]].xp[0]+node[indx[9]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.125*(node[indx[11]].xp[1]+node[indx[14]].xp[1]+node[indx[17]].xp[1]+node[indx[16]].xp[1]+
		      node[indx[15]].xp[1]+node[indx[12]].xp[1]+node[indx[9]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.125*(node[indx[11]].xp[2]+node[indx[14]].xp[2]+node[indx[17]].xp[2]+node[indx[16]].xp[2]+
		      node[indx[15]].xp[2]+node[indx[12]].xp[2]+node[indx[9]].xp[2]+node[indx[10]].xp[2]);
    }


    else if (node[v].ip==0 && node[v].ip!=gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp==0 && node[v].kp!=gd.LIMZ) { 
      node[v].x0[0]=0.5*(node[indx[16]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.5*(node[indx[16]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.5*(node[indx[16]].xp[2]+node[indx[10]].xp[2]);
    }
    else if (node[v].ip==0 && node[v].ip!=gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp!=0 && node[v].kp==gd.LIMZ) { 
      node[v].x0[0]=0.5*(node[indx[16]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.5*(node[indx[16]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.5*(node[indx[16]].xp[2]+node[indx[10]].xp[2]);
    }

    else if (node[v].ip!=0 && node[v].ip==gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp==0 && node[v].kp!=gd.LIMZ) { 
      node[v].x0[0]=0.5*(node[indx[16]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.5*(node[indx[16]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.5*(node[indx[16]].xp[2]+node[indx[10]].xp[2]);
    }
    else if (node[v].ip!=0 && node[v].ip==gd.LIMX && node[v].jp!=0 && node[v].jp!=gd.LIMY && node[v].kp!=0 && node[v].kp==gd.LIMZ) { 
      node[v].x0[0]=0.5*(node[indx[16]].xp[0]+node[indx[10]].xp[0]);
      node[v].x0[1]=0.5*(node[indx[16]].xp[1]+node[indx[10]].xp[1]);
      node[v].x0[2]=0.5*(node[indx[16]].xp[2]+node[indx[10]].xp[2]);
    }



    //loop over any Lagrangian interface J-Line specifications
    for (j=0;j<rlx.nlag;j++) {
      if (node[v].jp==rlx.jlines[j]) {
	node[v].x0[0]=node[v].xp[0];
	node[v].x0[1]=node[v].xp[1];
	node[v].x0[2]=node[v].xp[2];
      }
    }


  } /* end vertex loop */



  /* loop over all verticies to update xp */
  for ( auto v : mesh.vertices() ) {

      for (j=0;j<3;j++) {
	node[v].xp[j]=node[v].x0[j];
      }

  }

  
}
