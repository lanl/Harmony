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
void complete_stencil( 
  mymesh_t & m
) {
  int flag=0;
  int j;
  int c1,c2,c3,c4;

  flag=0;
  for (j=0;j<27;j++) 
    if (m.indx[j]!=-1) {
      m.cpy[j]=m.indx[j];
    }
    else
      flag=1;
	
  // fill in the coordinates (xc,yc,zc) corresponding to the -1 entries in the indx array
  // these are the ghost vertices
  if (flag==1) { //flag == 1 indicates boundary vertex      
    //cell cooresponding to +y face
    if (m.indx[16]==-1) {
      c2 = 13;
      c1 = 10;
      m.xc[16]=2.*m.xc[c2]-m.xc[c1];
      m.yc[16]=2.*m.yc[c2]-m.yc[c1];
      m.zc[16]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[16]=m.indx[13];
    }

    //cell cooresponding to -y face
    if (m.indx[10]==-1) {
      c2 = 13;
      c1 = 16;
      m.xc[10]=2.*m.xc[c2]-m.xc[c1];
      m.yc[10]=2.*m.yc[c2]-m.yc[c1];
      m.zc[10]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[10]=m.indx[13];
    }

    //cell cooresponding to +x face
    if (m.indx[14]==-1) {
      c2 = 13;
      c1 = 12;
      m.xc[14]=2.*m.xc[c2]-m.xc[c1];
      m.yc[14]=2.*m.yc[c2]-m.yc[c1];
      m.zc[14]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[14]=m.indx[13];
    }

    //cell cooresponding to -x face
    if (m.indx[12]==-1) {
      c2 = 13;
      c1 = 14;
      m.xc[12]=2.*m.xc[c2]-m.xc[c1];
      m.yc[12]=2.*m.yc[c2]-m.yc[c1];
      m.zc[12]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[12]=m.indx[13];
    }

    //cell cooresponding to the +z face
    if (m.indx[22]==-1) {
      c2 = 13;
      c1 = 4;
      m.xc[22]=2.*m.xc[c2]-m.xc[c1];
      m.yc[22]=2.*m.yc[c2]-m.yc[c1];
      m.zc[22]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[22]=m.indx[13];
    }

    //cell cooresponding to the -z face
    if (m.indx[4]==-1) {
      c2 = 13;
      c1 = 22;
      m.xc[4]=2.*m.xc[c2]-m.xc[c1];
      m.yc[4]=2.*m.yc[c2]-m.yc[c1];
      m.zc[4]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[4]=m.indx[13];
    }

    //the case when two face neighbors of the central cell in (z=0) plane are both boundary cells
    //top right cell
    if ((m.indx[16]==-1) && (m.indx[14]==-1)) {
      c3 = 12;
      c2 = 13;
      c1 = 10;
      m.xc[17]=2.*m.xc[c2]-m.xc[c1];
      m.yc[17]=2.*m.yc[c2]-m.yc[c1];
      m.zc[17]=2.*m.zc[c2]-m.zc[c1];
      m.xc[17]+=(m.xc[c2]-m.xc[c3]);
      m.yc[17]+=(m.yc[c2]-m.yc[c3]);
      m.zc[17]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[17]=m.indx[13];
    }
    else if ((m.indx[16]==-1) && (m.indx[17]==-1)) {
      c2 = 14;
      c1 = 11;
      m.xc[17]=2.*m.xc[c2]-m.xc[c1];
      m.yc[17]=2.*m.yc[c2]-m.yc[c1];
      m.zc[17]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[17]=m.indx[14];
    }
    else if ((m.indx[14]==-1) && (m.indx[17]==-1)) {
      c2 = 16;
      c1 = 15;
      m.xc[17]=2.*m.xc[c2]-m.xc[c1];
      m.yc[17]=2.*m.yc[c2]-m.yc[c1];
      m.zc[17]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[17]=m.indx[16];
    }

    //top left cell
    if ((m.indx[16]==-1) && (m.indx[12]==-1)) {
      c3 = 14;
      c2 = 13;
      c1 = 10;
      m.xc[15]=2.*m.xc[c2]-m.xc[c1];
      m.yc[15]=2.*m.yc[c2]-m.yc[c1];
      m.zc[15]=2.*m.zc[c2]-m.zc[c1];
      m.xc[15]+=(m.xc[c2]-m.xc[c3]);
      m.yc[15]+=(m.yc[c2]-m.yc[c3]);
      m.zc[15]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[15]=m.indx[13];
    }
    else if ((m.indx[16]==-1) && (m.indx[15]==-1)) {
      c2 = 12;
      c1 = 9;
      m.xc[15]=2.*m.xc[c2]-m.xc[c1];
      m.yc[15]=2.*m.yc[c2]-m.yc[c1];
      m.zc[15]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[15]=m.indx[12];
    }
    else if ((m.indx[12]==-1) && (m.indx[15]==-1)) {
      c2 = 16;
      c1 = 17;
      m.xc[15]=2.*m.xc[c2]-m.xc[c1];
      m.yc[15]=2.*m.yc[c2]-m.yc[c1];
      m.zc[15]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[15]=m.indx[16];
    }

    //bottom left cell
    if ((m.indx[10]==-1) && (m.indx[12]==-1)) {
      c3 = 14;
      c2 = 13;
      c1 = 16;
      m.xc[9]=2.*m.xc[c2]-m.xc[c1];
      m.yc[9]=2.*m.yc[c2]-m.yc[c1];
      m.zc[9]=2.*m.zc[c2]-m.zc[c1];
      m.xc[9]+=(m.xc[c2]-m.xc[c3]);
      m.yc[9]+=(m.yc[c2]-m.yc[c3]);
      m.zc[9]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[9]=m.indx[13];
    }
    else if ((m.indx[10]==-1) && (m.indx[9]==-1)) {
      c2 = 12;
      c1 = 15;
      m.xc[9]=2.*m.xc[c2]-m.xc[c1];
      m.yc[9]=2.*m.yc[c2]-m.yc[c1];
      m.zc[9]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[9]=m.indx[12];
    }
    else if ((m.indx[12]==-1) && (m.indx[9]==-1)) {
      c2 = 10;
      c1 = 11;
      m.xc[9]=2.*m.xc[c2]-m.xc[c1];
      m.yc[9]=2.*m.yc[c2]-m.yc[c1];
      m.zc[9]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[9]=m.indx[10];
    }

    //bottom right cell
    if ((m.indx[10]==-1) && (m.indx[14]==-1)) {
      c3 = 12;
      c2 = 13;
      c1 = 16;
      m.xc[11]=2.*m.xc[c2]-m.xc[c1];
      m.yc[11]=2.*m.yc[c2]-m.yc[c1];
      m.zc[11]=2.*m.zc[c2]-m.zc[c1];
      m.xc[11]+=(m.xc[c2]-m.xc[c3]);
      m.yc[11]+=(m.yc[c2]-m.yc[c3]);
      m.zc[11]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[11]=m.indx[13];
    }
    else if ((m.indx[10]==-1) && (m.indx[11]==-1)) {
      c2 = 14;
      c1 = 17;
      m.xc[11]=2.*m.xc[c2]-m.xc[c1];
      m.yc[11]=2.*m.yc[c2]-m.yc[c1];
      m.zc[11]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[11]=m.indx[14];
    }
    else if ((m.indx[14]==-1) && (m.indx[11]==-1)) {
      c2 = 10;
      c1 = 9;
      m.xc[11]=2.*m.xc[c2]-m.xc[c1];
      m.yc[11]=2.*m.yc[c2]-m.yc[c1];
      m.zc[11]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[11]=m.indx[10];
    }
    //end of the case when two face neighbors of the central cell in (z=0) plane are both boundary cells

    // +z plane
    //top right corner
    if ((m.indx[22]==-1) && (m.indx[16]==-1) && (m.indx[14]==-1)) {
      c4 = 4;
      c3 = 10;
      c2 = 13;
      c1 = 12;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.xc[26]+=(m.xc[c2]-m.xc[c3]);
      m.yc[26]+=(m.yc[c2]-m.yc[c3]);
      m.zc[26]+=(m.zc[c2]-m.zc[c3]);
      m.xc[26]+=(m.xc[c2]-m.xc[c4]);
      m.yc[26]+=(m.yc[c2]-m.yc[c4]);
      m.zc[26]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[26]=m.indx[13];
    }
    else if ((m.indx[16]!=-1) && (m.indx[25]==-1) && (m.indx[17]==-1)) {
      c3 = 7;
      c2 = 16;
      c1 = 15;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.xc[26]+=(m.xc[c2]-m.xc[c3]);
      m.yc[26]+=(m.yc[c2]-m.yc[c3]);
      m.zc[26]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[26]=m.indx[16];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[25]==-1) && (m.indx[23]==-1)) {
      c3 = 19;
      c2 = 22;
      c1 = 21;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.xc[26]+=(m.xc[c2]-m.xc[c3]);
      m.yc[26]+=(m.yc[c2]-m.yc[c3]);
      m.zc[26]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[26]=m.indx[22];
    }
    else if ((m.indx[14]!=-1) && (m.indx[17]==-1) && (m.indx[23]==-1)) {
      c3 = 5;
      c2 = 14;
      c1 = 11;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.xc[26]+=(m.xc[c2]-m.xc[c3]);
      m.yc[26]+=(m.yc[c2]-m.yc[c3]);
      m.zc[26]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[26]=m.indx[14];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[25]==-1) && (m.indx[26]==-1)) {
      c2 = 23;
      c1 = 20;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[26]=m.indx[23];
    }
    else if ((m.indx[22]!=-1) && (m.indx[23]==-1) && (m.indx[26]==-1)) {
      c2 = 25;
      c1 = 24;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[26]=m.indx[25];
    }
    else if ((m.indx[22]==-1) && (m.indx[26]==-1)) {
      c2 = 17;
      c1 = 8;
      m.xc[26]=2.*m.xc[c2]-m.xc[c1];
      m.yc[26]=2.*m.yc[c2]-m.yc[c1];
      m.zc[26]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[26]=m.indx[17];
    }

    //top 
    if ((m.indx[22]==-1) && (m.indx[25]==-1) && (m.indx[16]==-1)) {
      c3 = 4;
      c2 = 13;
      c1 = 10;
      m.xc[25]=2.*m.xc[c2]-m.xc[c1];
      m.yc[25]=2.*m.yc[c2]-m.yc[c1];
      m.zc[25]=2.*m.zc[c2]-m.zc[c1];
      m.xc[25]+=(m.xc[c2]-m.xc[c3]);
      m.yc[25]+=(m.yc[c2]-m.yc[c3]);
      m.zc[25]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[25]=m.indx[13];
    }
    else if ((m.indx[22]!=-1) && (m.indx[25]==-1) && (m.indx[26]==-1)) {
      c2 = 22;
      c1 = 19;
      m.xc[25]=2.*m.xc[c2]-m.xc[c1];
      m.yc[25]=2.*m.yc[c2]-m.yc[c1];
      m.zc[25]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[25]=m.indx[22];
    }
    else if ((m.indx[22]==-1) && (m.indx[25]==-1)) {
      c2 = 16;
      c1 = 7;
      m.xc[25]=2.*m.xc[c2]-m.xc[c1];
      m.yc[25]=2.*m.yc[c2]-m.yc[c1];
      m.zc[25]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[25]=m.indx[16];
    }

    //top left corner
    if ((m.indx[22]==-1) && (m.indx[16]==-1) && (m.indx[12]==-1)) {
      c4 = 4;
      c3 = 10;
      c2 = 13;
      c1 = 14;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.xc[24]+=(m.xc[c2]-m.xc[c3]);
      m.yc[24]+=(m.yc[c2]-m.yc[c3]);
      m.zc[24]+=(m.zc[c2]-m.zc[c3]);
      m.xc[24]+=(m.xc[c2]-m.xc[c4]);
      m.yc[24]+=(m.yc[c2]-m.yc[c4]);
      m.zc[24]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[24]=m.indx[13];
    }
    else if ((m.indx[16]!=-1) && (m.indx[25]==-1) && (m.indx[15]==-1)) {
      c3 = 7;
      c2 = 16;
      c1 = 17;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.xc[24]+=(m.xc[c2]-m.xc[c3]);
      m.yc[24]+=(m.yc[c2]-m.yc[c3]);
      m.zc[24]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[24]=m.indx[16];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[25]==-1) && (m.indx[21]==-1)) {
      c3 = 19;
      c2 = 22;
      c1 = 23;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.xc[24]+=(m.xc[c2]-m.xc[c3]);
      m.yc[24]+=(m.yc[c2]-m.yc[c3]);
      m.zc[24]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[24]=m.indx[22];
    }
    else if ((m.indx[12]!=-1) && (m.indx[21]==-1) && (m.indx[15]==-1)) {
      c3 = 3;
      c2 = 12;
      c1 = 9;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.xc[24]+=(m.xc[c2]-m.xc[c3]);
      m.yc[24]+=(m.yc[c2]-m.yc[c3]);
      m.zc[24]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[24]=m.indx[12];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[25]==-1) && (m.indx[24]==-1)) {
      c2 = 21;
      c1 = 18;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[24]=m.indx[21];
    }
    else if ((m.indx[22]!=-1) && (m.indx[21]==-1) && (m.indx[24]==-1)) {
      c2 = 25;
      c1 = 26;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[24]=m.indx[25];
    }
    else if ((m.indx[22]==-1) && (m.indx[24]==-1)) {
      c2 = 15;
      c1 = 6;
      m.xc[24]=2.*m.xc[c2]-m.xc[c1];
      m.yc[24]=2.*m.yc[c2]-m.yc[c1];
      m.zc[24]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[24]=m.indx[15];
    }

    //left 
    if ((m.indx[22]==-1) && (m.indx[21]==-1) && (m.indx[12]==-1)) {
      c3 = 4;
      c2 = 13;
      c1 = 14;
      m.xc[21]=2.*m.xc[c2]-m.xc[c1];
      m.yc[21]=2.*m.yc[c2]-m.yc[c1];
      m.zc[21]=2.*m.zc[c2]-m.zc[c1];
      m.xc[21]+=(m.xc[c2]-m.xc[c3]);
      m.yc[21]+=(m.yc[c2]-m.yc[c3]);
      m.zc[21]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[21]=m.indx[13];
    }
    else if ((m.indx[22]!=-1) && (m.indx[21]==-1) && (m.indx[24]==-1)) {
      c2 = 22;
      c1 = 23;
      m.xc[21]=2.*m.xc[c2]-m.xc[c1];
      m.yc[21]=2.*m.yc[c2]-m.yc[c1];
      m.zc[21]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[21]=m.indx[22];
    }
    else if ((m.indx[22]==-1) && (m.indx[21]==-1)) {
      c2 = 12;
      c1 = 3;
      m.xc[21]=2.*m.xc[c2]-m.xc[c1];
      m.yc[21]=2.*m.yc[c2]-m.yc[c1];
      m.zc[21]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[21]=m.indx[12];
    }

    //bottom left corner
    if ((m.indx[22]==-1) && (m.indx[10]==-1) && (m.indx[12]==-1)) {
      c4 = 4;
      c3 = 16;
      c2 = 13;
      c1 = 14;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.xc[18]+=(m.xc[c2]-m.xc[c3]);
      m.yc[18]+=(m.yc[c2]-m.yc[c3]);
      m.zc[18]+=(m.zc[c2]-m.zc[c3]);
      m.xc[18]+=(m.xc[c2]-m.xc[c4]);
      m.yc[18]+=(m.yc[c2]-m.yc[c4]);
      m.zc[18]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[18]=m.indx[13];
    }
    else if ((m.indx[10]!=-1) && (m.indx[19]==-1) && (m.indx[9]==-1)) {
      c3 = 1;
      c2 = 10;
      c1 = 11;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.xc[18]+=(m.xc[c2]-m.xc[c3]);
      m.yc[18]+=(m.yc[c2]-m.yc[c3]);
      m.zc[18]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[18]=m.indx[10];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[19]==-1) && (m.indx[21]==-1)) {
      c3 = 25;
      c2 = 22;
      c1 = 23;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.xc[18]+=(m.xc[c2]-m.xc[c3]);
      m.yc[18]+=(m.yc[c2]-m.yc[c3]);
      m.zc[18]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[18]=m.indx[22];
    }
    else if ((m.indx[12]!=-1) && (m.indx[9]==-1) && (m.indx[21]==-1)) {
      c3 = 3;
      c2 = 12;
      c1 = 15;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.xc[18]+=(m.xc[c2]-m.xc[c3]);
      m.yc[18]+=(m.yc[c2]-m.yc[c3]);
      m.zc[18]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[18]=m.indx[12];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[19]==-1) && (m.indx[18]==-1)) {
      c2 = 21;
      c1 = 24;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[18]=m.indx[21];
    }
    else if ((m.indx[22]!=-1) && (m.indx[21]==-1) && (m.indx[18]==-1)) {
      c2 = 19;
      c1 = 20;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[18]=m.indx[19];
    }
    else if ((m.indx[22]==-1) && (m.indx[18]==-1)) {
      c2 = 9;
      c1 = 0;
      m.xc[18]=2.*m.xc[c2]-m.xc[c1];
      m.yc[18]=2.*m.yc[c2]-m.yc[c1];
      m.zc[18]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[18]=m.indx[9];
    }

    //bottom
    if ((m.indx[22]==-1) && (m.indx[19]==-1) && (m.indx[10]==-1)) {
      c3 = 4;
      c2 = 13;
      c1 = 16;
      m.xc[19]=2.*m.xc[c2]-m.xc[c1];
      m.yc[19]=2.*m.yc[c2]-m.yc[c1];
      m.zc[19]=2.*m.zc[c2]-m.zc[c1];
      m.xc[19]+=(m.xc[c2]-m.xc[c3]);
      m.yc[19]+=(m.yc[c2]-m.yc[c3]);
      m.zc[19]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[19]=m.indx[13];
    }
    else if ((m.indx[22]!=-1) && (m.indx[19]==-1) && (m.indx[18]==-1)) {
      c2 = 22;
      c1 = 25;
      m.xc[19]=2.*m.xc[c2]-m.xc[c1];
      m.yc[19]=2.*m.yc[c2]-m.yc[c1];
      m.zc[19]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[19]=m.indx[22];
    }
    else if ((m.indx[22]==-1) && (m.indx[19]==-1)) {
      c2 = 10;
      c1 = 1;
      m.xc[19]=2.*m.xc[c2]-m.xc[c1];
      m.yc[19]=2.*m.yc[c2]-m.yc[c1];
      m.zc[19]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[19]=m.indx[10];
    }

    //bottom right corner
    if ((m.indx[22]==-1) && (m.indx[10]==-1) && (m.indx[14]==-1)) {
      c4 = 4;
      c3 = 16;
      c2 = 13;
      c1 = 12;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.xc[20]+=(m.xc[c2]-m.xc[c3]);
      m.yc[20]+=(m.yc[c2]-m.yc[c3]);
      m.zc[20]+=(m.zc[c2]-m.zc[c3]);
      m.xc[20]+=(m.xc[c2]-m.xc[c4]);
      m.yc[20]+=(m.yc[c2]-m.yc[c4]);
      m.zc[20]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[20]=m.indx[13];
    }
    else if ((m.indx[10]!=-1) && (m.indx[19]==-1) && (m.indx[11]==-1)) {
      c3 = 1;
      c2 = 10;
      c1 = 9;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.xc[20]+=(m.xc[c2]-m.xc[c3]);
      m.yc[20]+=(m.yc[c2]-m.yc[c3]);
      m.zc[20]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[20]=m.indx[10];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[19]==-1) && (m.indx[23]==-1)) {
      c3 = 25;
      c2 = 22;
      c1 = 21;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.xc[20]+=(m.xc[c2]-m.xc[c3]);
      m.yc[20]+=(m.yc[c2]-m.yc[c3]);
      m.zc[20]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[20]=m.indx[22];
    }
    else if ((m.indx[14]!=-1) && (m.indx[11]==-1) && (m.indx[23]==-1)) {
      c3 = 5;
      c2 = 14;
      c1 = 17;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.xc[20]+=(m.xc[c2]-m.xc[c3]);
      m.yc[20]+=(m.yc[c2]-m.yc[c3]);
      m.zc[20]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[20]=m.indx[14];
    }	
    else if ((m.indx[22]!=-1) && (m.indx[19]==-1) && (m.indx[20]==-1)) {
      c2 = 23;
      c1 = 26;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[20]=m.indx[23];
    }
    else if ((m.indx[22]!=-1) && (m.indx[23]==-1) && (m.indx[20]==-1)) {
      c2 = 19;
      c1 = 18;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[20]=m.indx[19];
    }
    else if ((m.indx[22]==-1) && (m.indx[20]==-1)) {
      c2 = 11;
      c1 = 2;
      m.xc[20]=2.*m.xc[c2]-m.xc[c1];
      m.yc[20]=2.*m.yc[c2]-m.yc[c1];
      m.zc[20]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[20]=m.indx[11];
    }

    //right
    if ((m.indx[22]==-1) && (m.indx[23]==-1) && (m.indx[14]==-1)) {
      c3 = 4;
      c2 = 13;
      c1 = 12;
      m.xc[23]=2.*m.xc[c2]-m.xc[c1];
      m.yc[23]=2.*m.yc[c2]-m.yc[c1];
      m.zc[23]=2.*m.zc[c2]-m.zc[c1];
      m.xc[23]+=(m.xc[c2]-m.xc[c3]);
      m.yc[23]+=(m.yc[c2]-m.yc[c3]);
      m.zc[23]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[23]=m.indx[13];
    }
    else if ((m.indx[22]!=-1) && (m.indx[23]==-1) && (m.indx[20]==-1)) {
      c2 = 22;
      c1 = 21;
      m.xc[23]=2.*m.xc[c2]-m.xc[c1];
      m.yc[23]=2.*m.yc[c2]-m.yc[c1];
      m.zc[23]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[23]=m.indx[22];
    }
    else if ((m.indx[22]==-1) && (m.indx[23]==-1)) {
      c2 = 14;
      c1 = 5;
      m.xc[23]=2.*m.xc[c2]-m.xc[c1];
      m.yc[23]=2.*m.yc[c2]-m.yc[c1];
      m.zc[23]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[23]=m.indx[14];
    }
    //end of +z plane

    // -z plane
    //top right corner
    if ((m.indx[4]==-1) && (m.indx[16]==-1) && (m.indx[14]==-1)) {
      c4 = 22;
      c3 = 10;
      c2 = 13;
      c1 = 12;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.xc[8]+=(m.xc[c2]-m.xc[c3]);
      m.yc[8]+=(m.yc[c2]-m.yc[c3]);
      m.zc[8]+=(m.zc[c2]-m.zc[c3]);
      m.xc[8]+=(m.xc[c2]-m.xc[c4]);
      m.yc[8]+=(m.yc[c2]-m.yc[c4]);
      m.zc[8]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[8]=m.indx[13];
    }
    else if ((m.indx[16]!=-1) && (m.indx[7]==-1) && (m.indx[17]==-1)) {
      c3 = 25;
      c2 = 16;
      c1 = 15;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.xc[8]+=(m.xc[c2]-m.xc[c3]);
      m.yc[8]+=(m.yc[c2]-m.yc[c3]);
      m.zc[8]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[8]=m.indx[16];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[7]==-1) && (m.indx[5]==-1)) {
      c3 = 1;
      c2 = 4;
      c1 = 3;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.xc[8]+=(m.xc[c2]-m.xc[c3]);
      m.yc[8]+=(m.yc[c2]-m.yc[c3]);
      m.zc[8]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[8]=m.indx[4];
    }
    else if ((m.indx[14]!=-1) && (m.indx[17]==-1) && (m.indx[5]==-1)) {
      c3 = 23;
      c2 = 14;
      c1 = 11;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.xc[8]+=(m.xc[c2]-m.xc[c3]);
      m.yc[8]+=(m.yc[c2]-m.yc[c3]);
      m.zc[8]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[8]=m.indx[14];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[7]==-1) && (m.indx[8]==-1)) {
      c2 = 5;
      c1 = 2;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[8]=m.indx[5];
    }
    else if ((m.indx[4]!=-1) && (m.indx[5]==-1) && (m.indx[8]==-1)) {
      c2 = 7;
      c1 = 6;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[8]=m.indx[7];
    }
    else if ((m.indx[4]==-1) && (m.indx[8]==-1)) {
      c2 = 17;
      c1 = 26;
      m.xc[8]=2.*m.xc[c2]-m.xc[c1];
      m.yc[8]=2.*m.yc[c2]-m.yc[c1];
      m.zc[8]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[8]=m.indx[17];
    }

    //top
    if ((m.indx[4]==-1) && (m.indx[7]==-1) && (m.indx[16]==-1)) {
      c3 = 22;
      c2 = 13;
      c1 = 10;
      m.xc[7]=2.*m.xc[c2]-m.xc[c1];
      m.yc[7]=2.*m.yc[c2]-m.yc[c1];
      m.zc[7]=2.*m.zc[c2]-m.zc[c1];
      m.xc[7]+=(m.xc[c2]-m.xc[c3]);
      m.yc[7]+=(m.yc[c2]-m.yc[c3]);
      m.zc[7]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[7]=m.indx[13];
    }
    else if ((m.indx[4]!=-1) && (m.indx[7]==-1) && (m.indx[8]==-1)) {
      c2 = 4;
      c1 = 1;
      m.xc[7]=2.*m.xc[c2]-m.xc[c1];
      m.yc[7]=2.*m.yc[c2]-m.yc[c1];
      m.zc[7]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[7]=m.indx[4];
    }
    else if ((m.indx[4]==-1) && (m.indx[7]==-1)) {
      c2 = 16;
      c1 = 25;
      m.xc[7]=2.*m.xc[c2]-m.xc[c1];
      m.yc[7]=2.*m.yc[c2]-m.yc[c1];
      m.zc[7]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[7]=m.indx[16];
    }

    //top left corner
    if ((m.indx[4]==-1) && (m.indx[16]==-1) && (m.indx[12]==-1)) {
      c4 = 22;
      c3 = 10;
      c2 = 13;
      c1 = 14;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.xc[6]+=(m.xc[c2]-m.xc[c3]);
      m.yc[6]+=(m.yc[c2]-m.yc[c3]);
      m.zc[6]+=(m.zc[c2]-m.zc[c3]);
      m.xc[6]+=(m.xc[c2]-m.xc[c4]);
      m.yc[6]+=(m.yc[c2]-m.yc[c4]);
      m.zc[6]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[6]=m.indx[13];
    }
    else if ((m.indx[16]!=-1) && (m.indx[7]==-1) && (m.indx[15]==-1)) {
      c3 = 25;
      c2 = 16;
      c1 = 17;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.xc[6]+=(m.xc[c2]-m.xc[c3]);
      m.yc[6]+=(m.yc[c2]-m.yc[c3]);
      m.zc[6]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[6]=m.indx[16];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[7]==-1) && (m.indx[3]==-1)) {
      c3 = 1;
      c2 = 4;
      c1 = 5;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.xc[6]+=(m.xc[c2]-m.xc[c3]);
      m.yc[6]+=(m.yc[c2]-m.yc[c3]);
      m.zc[6]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[6]=m.indx[4];
    }
    else if ((m.indx[12]!=-1) && (m.indx[15]==-1) && (m.indx[3]==-1)) {
      c3 = 21;
      c2 = 12;
      c1 = 9;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.xc[6]+=(m.xc[c2]-m.xc[c3]);
      m.yc[6]+=(m.yc[c2]-m.yc[c3]);
      m.zc[6]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[6]=m.indx[12];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[7]==-1) && (m.indx[6]==-1)) {
      c2 = 3;
      c1 = 0;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[6]=m.indx[3];
    }
    else if ((m.indx[4]!=-1) && (m.indx[3]==-1) && (m.indx[6]==-1)) {
      c2 = 7;
      c1 = 8;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[6]=m.indx[7];
    }
    else if ((m.indx[4]==-1) && (m.indx[6]==-1)) {
      c2 = 15;
      c1 = 24;
      m.xc[6]=2.*m.xc[c2]-m.xc[c1];
      m.yc[6]=2.*m.yc[c2]-m.yc[c1];
      m.zc[6]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[6]=m.indx[15];
    }

    //left 
    if ((m.indx[4]==-1) && (m.indx[3]==-1) && (m.indx[12]==-1)) {
      c3 = 22;
      c2 = 13;
      c1 = 14;
      m.xc[3]=2.*m.xc[c2]-m.xc[c1];
      m.yc[3]=2.*m.yc[c2]-m.yc[c1];
      m.zc[3]=2.*m.zc[c2]-m.zc[c1];
      m.xc[3]+=(m.xc[c2]-m.xc[c3]);
      m.yc[3]+=(m.yc[c2]-m.yc[c3]);
      m.zc[3]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[3]=m.indx[13];
    }
    else if ((m.indx[4]!=-1) && (m.indx[3]==-1) && (m.indx[6]==-1)) {
      c2 = 4;
      c1 = 5;
      m.xc[3]=2.*m.xc[c2]-m.xc[c1];
      m.yc[3]=2.*m.yc[c2]-m.yc[c1];
      m.zc[3]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[3]=m.indx[4];
    }
    else if ((m.indx[4]==-1) && (m.indx[3]==-1)) {
      c2 = 12;
      c1 = 21;
      m.xc[3]=2.*m.xc[c2]-m.xc[c1];
      m.yc[3]=2.*m.yc[c2]-m.yc[c1];
      m.zc[3]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[3]=m.indx[12];
    }
    
    //bottom left corner
    if ((m.indx[4]==-1) && (m.indx[10]==-1) && (m.indx[12]==-1)) {
      c4 = 22;
      c3 = 16;
      c2 = 13;
      c1 = 14;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.xc[0]+=(m.xc[c2]-m.xc[c3]);
      m.yc[0]+=(m.yc[c2]-m.yc[c3]);
      m.zc[0]+=(m.zc[c2]-m.zc[c3]);
      m.xc[0]+=(m.xc[c2]-m.xc[c4]);
      m.yc[0]+=(m.yc[c2]-m.yc[c4]);
      m.zc[0]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[0]=m.indx[13];
    }
    else if ((m.indx[10]!=-1) && (m.indx[1]==-1) && (m.indx[9]==-1)) {
      c3 = 19;
      c2 = 10;
      c1 = 11;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.xc[0]+=(m.xc[c2]-m.xc[c3]);
      m.yc[0]+=(m.yc[c2]-m.yc[c3]);
      m.zc[0]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[0]=m.indx[10];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[1]==-1) && (m.indx[3]==-1)) {
      c3 = 7;
      c2 = 4;
      c1 = 5;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.xc[0]+=(m.xc[c2]-m.xc[c3]);
      m.yc[0]+=(m.yc[c2]-m.yc[c3]);
      m.zc[0]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[0]=m.indx[4];
    }
    else if ((m.indx[12]!=-1) && (m.indx[9]==-1) && (m.indx[3]==-1)) {
      c3 = 21;
      c2 = 12;
      c1 = 15;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.xc[0]+=(m.xc[c2]-m.xc[c3]);
      m.yc[0]+=(m.yc[c2]-m.yc[c3]);
      m.zc[0]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[0]=m.indx[12];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[1]==-1) && (m.indx[0]==-1)) {
      c2 = 3;
      c1 = 6;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[0]=m.indx[3];
    }
    else if ((m.indx[4]!=-1) && (m.indx[3]==-1) && (m.indx[0]==-1)) {
      c2 = 1;
      c1 = 2;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[0]=m.indx[1];
    }
    else if ((m.indx[4]==-1) && (m.indx[0]==-1)) {
      c2 = 9;
      c1 = 18;
      m.xc[0]=2.*m.xc[c2]-m.xc[c1];
      m.yc[0]=2.*m.yc[c2]-m.yc[c1];
      m.zc[0]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[0]=m.indx[9];
    }

    //bottom
    if ((m.indx[4]==-1) && (m.indx[1]==-1) && (m.indx[10]==-1)) {
      c3 = 22;
      c2 = 13;
      c1 = 16;
      m.xc[1]=2.*m.xc[c2]-m.xc[c1];
      m.yc[1]=2.*m.yc[c2]-m.yc[c1];
      m.zc[1]=2.*m.zc[c2]-m.zc[c1];
      m.xc[1]+=(m.xc[c2]-m.xc[c3]);
      m.yc[1]+=(m.yc[c2]-m.yc[c3]);
      m.zc[1]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[1]=m.indx[13];
    }
    else if ((m.indx[4]!=-1) && (m.indx[1]==-1) && (m.indx[0]==-1)) {
      c2 = 4;
      c1 = 7;
      m.xc[1]=2.*m.xc[c2]-m.xc[c1];
      m.yc[1]=2.*m.yc[c2]-m.yc[c1];
      m.zc[1]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[1]=m.indx[4];
    }
    else if ((m.indx[4]==-1) && (m.indx[1]==-1)) {
      c2 = 10;
      c1 = 19;
      m.xc[1]=2.*m.xc[c2]-m.xc[c1];
      m.yc[1]=2.*m.yc[c2]-m.yc[c1];
      m.zc[1]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[1]=m.indx[10];
    }

    //bottom right corner
    if ((m.indx[4]==-1) && (m.indx[10]==-1) && (m.indx[14]==-1)) {
      c4 = 22;
      c3 = 16;
      c2 = 13;
      c1 = 12;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.xc[2]+=(m.xc[c2]-m.xc[c3]);
      m.yc[2]+=(m.yc[c2]-m.yc[c3]);
      m.zc[2]+=(m.zc[c2]-m.zc[c3]);
      m.xc[2]+=(m.xc[c2]-m.xc[c4]);
      m.yc[2]+=(m.yc[c2]-m.yc[c4]);
      m.zc[2]+=(m.zc[c2]-m.zc[c4]);
      m.cpy[2]=m.indx[13];
    }
    else if ((m.indx[10]!=-1) && (m.indx[1]==-1) && (m.indx[11]==-1)) {
      c3 = 19;
      c2 = 10;
      c1 = 9;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.xc[2]+=(m.xc[c2]-m.xc[c3]);
      m.yc[2]+=(m.yc[c2]-m.yc[c3]);
      m.zc[2]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[2]=m.indx[10];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[1]==-1) && (m.indx[5]==-1)) {
      c3 = 7;
      c2 = 4;
      c1 = 3;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.xc[2]+=(m.xc[c2]-m.xc[c3]);
      m.yc[2]+=(m.yc[c2]-m.yc[c3]);
      m.zc[2]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[2]=m.indx[4];
    }
    else if ((m.indx[14]!=-1) && (m.indx[11]==-1) && (m.indx[5]==-1)) {
      c3 = 23;
      c2 = 14;
      c1 = 17;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.xc[2]+=(m.xc[c2]-m.xc[c3]);
      m.yc[2]+=(m.yc[c2]-m.yc[c3]);
      m.zc[2]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[2]=m.indx[14];
    }	
    else if ((m.indx[4]!=-1) && (m.indx[1]==-1) && (m.indx[2]==-1)) {
      c2 = 5;
      c1 = 8;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[2]=m.indx[5];
    }
    else if ((m.indx[4]!=-1) && (m.indx[5]==-1) && (m.indx[2]==-1)) {
      c2 = 1;
      c1 = 0;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[2]=m.indx[1];
    }
    else if ((m.indx[4]==-1) && (m.indx[2]==-1)) {
      c2 = 11;
      c1 = 20;
      m.xc[2]=2.*m.xc[c2]-m.xc[c1];
      m.yc[2]=2.*m.yc[c2]-m.yc[c1];
      m.zc[2]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[2]=m.indx[11];
    }

    //right
    if ((m.indx[4]==-1) && (m.indx[5]==-1) && (m.indx[14]==-1)) {
      c3 = 22;
      c2 = 13;
      c1 = 12;
      m.xc[5]=2.*m.xc[c2]-m.xc[c1];
      m.yc[5]=2.*m.yc[c2]-m.yc[c1];
      m.zc[5]=2.*m.zc[c2]-m.zc[c1];
      m.xc[5]+=(m.xc[c2]-m.xc[c3]);
      m.yc[5]+=(m.yc[c2]-m.yc[c3]);
      m.zc[5]+=(m.zc[c2]-m.zc[c3]);
      m.cpy[5]=m.indx[13];
    }
    else if ((m.indx[4]!=-1) && (m.indx[5]==-1) && (m.indx[2]==-1)) {
      c2 = 4;
      c1 = 3;
      m.xc[5]=2.*m.xc[c2]-m.xc[c1];
      m.yc[5]=2.*m.yc[c2]-m.yc[c1];
      m.zc[5]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[5]=m.indx[4];
    }
    else if ((m.indx[4]==-1) && (m.indx[5]==-1)) {
      c2 = 14;
      c1 = 23;
      m.xc[5]=2.*m.xc[c2]-m.xc[c1];
      m.yc[5]=2.*m.yc[c2]-m.yc[c1];
      m.zc[5]=2.*m.zc[c2]-m.zc[c1];
      m.cpy[5]=m.indx[14];
    }
      
    //end of -z plane
  } //end if flag
   
}
