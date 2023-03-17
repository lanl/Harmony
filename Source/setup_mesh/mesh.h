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
#define MAX(a,b) ((b)>(a) ? (b) : (a))
#define MIN(a,b) ((b)>(a) ? (a) : (b))

int findex(int max, int i, int j) {

  return max*(i-1)+j-1;

}

int findex1(int i) {

  return (i-1);

}

class mesh_t {
  
public :
  
  int kk3ll;
  int igeo;
  int kk1ll,kk1a,kk1v;
  int kk6ll,k00,k11,k22,k12,k13,k23;
  int kk9ll,k01,k02,k03;

  int kkpl, kkpll, kkpln;
  int kkzl, kkzll, kkzln;
  int kksll, kksl, kksln;
  int kkfll, kkfl;
  int kkcll, kkell, kkcl, kkel;
  int kkgll;

  int *kptyp; 
  int *kstyp;
  int *kksz;
  int *kkss4;
  int *kksp1;
  int *kksp2;
  int *kkss2_ejl;


  double *px; //2D array
  int *indx; //2D array

 
  int *kkss5;
  int *kspe;
  int *kkss2_ejlpe;
  int *kkss2;
  int *kksz2_ejlpe;
  int *kksz2_ejl;
  int *kkss3;
  int *kksz2;
  int *kztyp;
  int *kkzs;
  int *ksleft;
  int *kkse;
  int *ketyp;
  int *kkep1;
  int *kkep2;
  int *kksf;
  int *kftyp;
  int *kfpe;
  int *kkfz2_ejl;
  int *kkfz2_ejlpe;
  int *kkfz2;
  int *kkfz;
  int *kksc1;
  int *kksc2;
  int *kkcz;
  int *kkcp;
  int *kctyp;
  int *kkcs1;
  int *kkcs2;

  mesh_t(char *buffer);

};


mesh_t::mesh_t(char *buffer)
{
  FILE *in;
  int i,j,k;
  double xa;
  int nz,ng,np;
  int maxping,maxginz;

  int g,g2,n,p,ip1,ip2,jp1,jp2,jp3,n2,ipass;
  int if0,if1,z,s,sn,nsides,m,pa,s2,maxgatp;
  int kkzlx,kkslx,num,z1,z2,p1,p2,pp1,pp2;

  int c,c1,c2,e,ss5,s0,max_znumf,f,s_now,s_next;

  char cc;

  int *ip; 
  int *ig;
  int *nping; 
  int *iping;
  int *igtyp;
  int *ising;
  int *nng2;
  int *nngz;
  int *iz;
  int *nginz;
  int *iginz;
  int *iztyp;
  int *ipinv;
  int *iginv;
  int *iptyp;  
  int *numgatp;
  int *listgatp;
  
  int *numz;
  int *nnzs;
  int *nump;
  int *nnpp;
  int *nnpe;
  int *kstag;
  int *nnzz;
  int *nnzzpe;
  int *nnzf;
  int *nnzp;
  int *nnzc;  

  
  kk3ll = 3;
  igeo=0;
 
  kk1ll=1;
  kk1v=1;
  kk1a=1;
 
  kk6ll=6;
  k00=1;
  k11=2;
  k22=3;
  k12=4;
  k13=5;
  k23=6;
 
  kk9ll=9;
  k01=7;
  k02=8;
  k03=9;

  in=fopen(buffer,"r");
  for (j=1;j<=12;j++) {

    i=0;
    while ((cc=(char)fgetc(in))!='\n') {
        i++;
    }
  }

  fscanf(in,"%d",&kkpll);

  px = new double[3*kkpll];
  kkpln=kkpll;
  
  for (i=1;i<=kkpll;i++) {

    fscanf(in,"%d",&k);

  }
    
  for (i=1;i<=kkpll;i++) {

    fscanf(in,"%le",&xa);
    px[findex(3,i,1)]=xa;

  }

  for (i=1;i<=kkpll;i++) {

    fscanf(in,"%le",&xa);
    px[findex(3,i,2)]=xa;

  }

  for (i=1;i<=kkpll;i++) {

    fscanf(in,"%le",&xa);
    px[findex(3,i,3)]=xa;

  }


  cc=(char)fgetc(in);
  for (j=1;j<=1;j++) {
    i=0;
    while ((cc=(char)fgetc(in))!='\n') {
      i++;
    }
  }

  fscanf(in,"%d",&kkzln);
  
  kkzll=kkzln;
  indx = new int[8*kkzll];

  

  for (i=1;i<=kkzln;i++) {
      
    fscanf(in,"%d",&k);

  } 
  
  for (i=1;i<=kkzln;i++) {

    j=findex(8,i,1);
    fscanf(in,"%d%d%d%d%d%d%d%d",indx+j,indx+j+1,indx+j+2,indx+j+3,
	   indx+j+4,indx+j+5,indx+j+6,indx+j+7);

  }

  fclose(in);
  
  maxping=4;
  maxginz=6;

  np = kkpll;
  nz = kkzll;
  ng = kkgll = nz*6;

  ip = new int[np];
  ig = new int[ng];
  nping = new int[ng];
  iping = new int[ng*maxping];
  igtyp = new int[ng];
  ising = new int[ng*maxping];
  nng2 = new int[ng];
  nngz = new int[ng];
  iz = new int[nz];
  nginz = new int[nz]; 
  iginz = new int[nz*maxginz];
  iztyp = new int[nz];  

  for (k=1;k<=np;k++) 
    ip[findex1(k)]=k;
  
  g = 1;

  for (k=1;k<=nz;k++) {

    //bottom face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,1)];
    iping[findex(4,g,2)]=indx[findex(8,k,4)];
    iping[findex(4,g,3)]=indx[findex(8,k,8)];
    iping[findex(4,g,4)]=indx[findex(8,k,5)];    
    igtyp[findex1(g)]=1;
    g++;
    
    //left face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,1)];
    iping[findex(4,g,2)]=indx[findex(8,k,5)];
    iping[findex(4,g,3)]=indx[findex(8,k,6)];
    iping[findex(4,g,4)]=indx[findex(8,k,2)];    
    igtyp[findex1(g)]=1;
    g++;    

    //front face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,1)];
    iping[findex(4,g,2)]=indx[findex(8,k,2)];
    iping[findex(4,g,3)]=indx[findex(8,k,3)];
    iping[findex(4,g,4)]=indx[findex(8,k,4)];    
    igtyp[findex1(g)]=1;
    g++;        
                    
    //right face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,4)];
    iping[findex(4,g,2)]=indx[findex(8,k,3)];
    iping[findex(4,g,3)]=indx[findex(8,k,7)];
    iping[findex(4,g,4)]=indx[findex(8,k,8)];    
    igtyp[findex1(g)]=1;
    g++;            
                  
    //back face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,6)];
    iping[findex(4,g,2)]=indx[findex(8,k,5)];
    iping[findex(4,g,3)]=indx[findex(8,k,8)];
    iping[findex(4,g,4)]=indx[findex(8,k,7)];    
    igtyp[findex1(g)]=1;
    g++;         
                
    //top face
    ig[findex1(g)]=g;
    nping[findex1(g)]=4;
    iping[findex(4,g,1)]=indx[findex(8,k,6)];
    iping[findex(4,g,2)]=indx[findex(8,k,7)];
    iping[findex(4,g,3)]=indx[findex(8,k,3)];
    iping[findex(4,g,4)]=indx[findex(8,k,2)];    
    igtyp[findex1(g)]=1;
    g++;       
  }


  for (k=1;k<=nz;k++) {
    iz[findex1(k)]=k;
    nginz[findex1(k)]=6;

    if0 = (k-1)*6 + 1;
    if1 = if0 + 5;
    n=1;
    for(j = if0; j<=if1;j++) {
      iginz[findex(6,k,n)]=j;
      n++;
    }
    iztyp[findex1(k)]=1;
  }
  
  ipinv= new int[kkpln];
  
  for (k=1;k<=np;k++) {  
    m=ip[findex1(k)];
    ipinv[findex1(m)]=ip[findex1(k)];
  }

  numgatp= new int[np];
  for (k=1;k<=np;k++)
    numgatp[findex1(k)]=0;

  maxgatp = 0;
  for (k=1;k<=ng;k++) {
    for (n=1;n<=nping[findex1(k)];n++) {
      p = iping[findex(4,k,n)];
      if (p==0) printf("g= %d n= %d p= %d\n",k,n,iping[findex(4,k,n)]);
      numgatp[findex1(p)] = numgatp[findex1(p)] + 1;
      maxgatp = MAX(maxgatp,numgatp[findex1(p)]);
    }
  }
  
  listgatp=new int[np*maxgatp];
  for (k=1;k<=np;k++) 
    for (n=1;n<=maxgatp;n++)
      listgatp[findex(maxgatp,k,n)]=0;

  for (g=1;g<=ng;g++) 
    for (n=1;n<=nping[findex1(g)];n++) {
      p = iping[findex(4,g,n)];
      for (i=1;i<=numgatp[findex1(p)];i++) {
	if (listgatp[findex(maxgatp,p,i)]==g)
	  break;
	else if (listgatp[findex(maxgatp,p,i)]==0){
	  listgatp[findex(maxgatp,p,i)]= g;
	  break;
	}
      }
    }
  
  for (g=1;g<=ng;g++) {
    nng2[findex1(g)] = 0;
    ip1     = iping[findex(4,g,1)];
    p1      = ipinv[findex1(ip1)];
    for (i=1;i<=numgatp[findex1(p1)];i++) {
      g2 = listgatp[findex(maxgatp,p1,i)];
      if (g2==0) printf("fatal error in face setup g= %d num= %d\n",g,numgatp[findex1(p1)]);
      if ((g2!=g) && (nping[findex1(g)]==nping[findex1(g2)])) {
	for (n=1;n<=nping[findex1(g2)];n++) {
	  k=(n % nping[findex1(g2)])+1;
	  jp1=iping[findex(4,g2,k)];
	  if (jp1==ip1) {
	    ipass=1;
	    for (j=1;j<=nping[findex1(g2)];j++) {
	      k=((n+j) % nping[findex1(g2)])+1;
	      jp2=iping[findex(4,g2,k)];
	      k=((nping[findex1(g)]-j) % nping[findex1(g)])+1;
	      ip2=iping[findex(4,g,k)];
	      if (jp2!=ip2) {
		ipass=0;
		break;
	      }
	    }
	    if (ipass==1) {
	      nng2[findex1(g)]=g2;
	      break;
	    }
	  }
	} //n
      }
    } //i
  }


  kksln = 0;
  kkzlx = 0;
  kkslx = 0;
  for (k=1;k<=ng;k++) {  
    nsides = nping[findex1(k)];
    kksln  = kksln + nsides;

    if (nng2[findex1(k)]==0) {
      kkslx = kkslx + nsides;
      kkzlx = kkzlx + 1;
    }
  }

  kksll = kksln + kkslx;
  kkzll = kkzln + kkzlx;

  kkpll = kkpln;

  kksl=kksll;
  kkzl=kkzll;
  kkpl=kkpll;


  //allocate side variables
  kptyp = new int[kkpll];
  for (k=1;k<=kkpll;k++)
    kptyp[findex1(k)]=0;
  kstyp = new int[kksll];
  kksz = new int[kksll];
  kkss4 = new int[kksll];
  kksp1 = new int[kksll];
  kksp2 = new int[kksll];
  kkss2_ejl = new int[kksll];
  for (k=1;k<=kksll;k++)
    kkss2_ejl[findex1(k)]=0;  

  iginv = new int[kkgll];
  iptyp = new int[np];


  for (k=1;k<=kkgll;k++)
    iginv[findex1(k)]=0;

  for (k=1;k<=np;k++)
    iptyp[findex1(k)]=1;

  for (g=1;g<=ng;g++) {
    m=ig[findex1(g)];
    iginv[findex1(m)]=g;
  }

  for (m=1;m<=nz;m++) 
    for (n=1;n<=nginz[findex1(m)];n++) {
      g = iginz[findex(6,m,n)];
      nngz[findex1(g)] = iz[findex1(m)];
    }

  s = 0;
  for (g=1;g<=ng;g++)
    for (n = 1;n<=nping[findex1(g)];n++) {
      s = s + 1;
      ising[findex(4,g,n)] = s;
      kstyp[findex1(s)] = 1;
      kksz[findex1(s)] = nngz[findex1(g)];
      kksp1[findex1(s)] = iping[findex(4,g,n)];
      k=(n % nping[findex1(g)])+1;
      kksp2[findex1(s)] = iping[findex(4,g,k)];
          
      if (n==nping[findex1(g)])
	kkss4[findex1(s)] = ising[findex(4,g,1)];
      else
	kkss4[findex1(s)] = s + 1;
    }

  z = kkzln;
  s = kksln;
  
  for (g=1;g<=ng;g++)
    if (nng2[findex1(g)]==0) {
      z  = z + 1;
      sn = s + 1;
      for (n=nping[findex1(g)];n>=1;n--) {
	s = s + 1;
	kstyp[findex1(s)]  = -1;
	kksz[findex1(s)]   = z;
	kksp2[findex1(s)]  = iping[findex(4,g,n)];
	kksp1[findex1(s)]  = iping[findex(4,g,(n % nping[findex1(g)])+1)];
	s2 = ising[findex(4,g,n)];
	kkss2_ejl[findex1(s)]  = s2;
	kkss2_ejl[findex1(s2)] = s;
	if (n==1) 
	  kkss4[findex1(s)] = sn;
	else
	  kkss4[findex1(s)] = s + 1;
      }
    }

  for (s=1;s<=kksl;s++) {
    p1=kksp1[findex1(s)];
    p2=kksp2[findex1(s)];
    if (kstyp[findex1(s)]<0) {
      kptyp[findex1(p1)]=-1;
      kptyp[findex1(p2)]=-1;
    }
    else if (kstyp[findex1(s)]>0) {
      if (kptyp[findex1(p1)]==0) kptyp[findex1(p1)]=1;
      if (kptyp[findex1(p2)]==0) kptyp[findex1(p2)]=1;
    }
  }

  for (k=1;k<=kkpl;k++)
    if (iptyp[findex1(k)]==0) kptyp[findex1(k)]=0;

  
  for (g=1;g<=ng;g++) {
    for (n=1;n<=nping[findex1(g)];n++) {
      s  = ising[findex(4,g,n)];
      pa = kksp1[findex1(s)];
      g2 = nng2[findex1(g)];
      if (g2!=0) {
	for (n2=1;n2<=nping[findex1(g2)];n2++) {
	  s2 = ising[findex(4,g2,n2)];
	  if (kksp2[findex1(s2)]==pa) {
	    if (kkss2_ejl[findex1(s)]==0) 
	      kkss2_ejl[findex1(s)] = s2;
	    else if (kkss2_ejl[findex1(s)]!=s2) 
              printf("kkss2_ejl error 1\n");

	    if (kkss2_ejl[findex1(s2)]==0) 
	      kkss2_ejl[findex1(s2)] = s;
	    else if (kkss2_ejl[findex1(s2)]!=s)
	      printf("kkss2_ejl error 2\n");

	    goto stop;
	  }
	} //n2
      }
					    
    stop: continue;
    } //n
  } //g
	
  delete[] ip; 
  delete[] ig;
  delete[] nping; 
  delete[] iping;
  delete[] igtyp;
  delete[] ising;
  delete[] nng2;
  delete[] nngz;
  delete[] iz;
  delete[] nginz;
  delete[] iginz;
  delete[] iztyp;
  delete[] ipinv;
  delete[] iginv;
  delete[] iptyp;  
  delete[] numgatp;
  delete[] listgatp;


  //setup sides
  kspe = new int[kksll];
  kkss2_ejlpe = new int[kksll];

  for (k=1;k<=kksll;k++) {
    kspe[findex1(k)]=0;
    kkss2_ejlpe[findex1(k)]=0;
  }

  kkss2= new int[kksll];
 
  for (s=1;s<=kksl;s++) 
    if (kstyp[findex1(s)]!=0) 
      kkss2[findex1(s)] = kkss2_ejl[findex1(s)];

  kksz2_ejlpe= new int[kksll];
  kksz2_ejl= new int[kksll];
  
  for (s=1;s<=kksl;s++)
    kksz2_ejlpe[findex1(s)] = kkss2_ejlpe[findex1(s)];

  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]!=0) {
      s2 = kkss2_ejl[findex1(s)];
      kksz2_ejl[findex1(s)] = kksz[findex1(s2)];
    }
  
  kkss3= new int[kksll];
  for (s0=1;s0<=kksl;s0++) { 
    s = s0;
    while (kkss4[findex1(s)]!=s0)
      s = kkss4[findex1(s)];

    kkss3[findex1(s0)] = s;
  }

  kksz2= new int[kksll];
  for (s=1;s<=kksll;s++)
    kksz2[findex1(s)]=0;
  
  for (s=1;s<=kksl;s++)    
    if(kstyp[findex1(s)]!=0)
      kksz2[findex1(s)] = kksz2_ejl[findex1(s)];
     
  kkzs= new int[kksll];
  for (s=1;s<=kksll;s++)
    kkzs[findex1(s)]=0;
  
  for (s=1;s<=kksl;s++)
    if(kstyp[findex1(s)]!=0) {
      z=kksz[findex1(s)];
      kkzs[findex1(z)]=s;
    }

  kztyp= new int[kkzll];
  for (z=1;z<=kkzll;z++)
    kztyp[findex1(z)]=0;
  
  for (z=1;z<=kkzl;z++) {
    s=kkzs[findex1(z)];
    if (s==0) 
      kztyp[findex1(z)]=0;
    else
      kztyp[findex1(z)]=MIN(1, kstyp[findex1(s)]);
  }

  numz= new int[kkzll];
  for (z=1;z<=kkzl;z++)
    numz[findex1(z)]=0;
  
  num=0;
  for (s=1;s<=kksl;s++)
    if(kstyp[findex1(s)]!=0)  { 
      z=kksz[findex1(s)];
      numz[findex1(z)]=numz[findex1(z)]+1;
      num=MAX(num,numz[findex1(z)]);
    }

  nnzs= new int[kkzll*num];
  for (z=1;z<=kkzll;z++)
    for (n=1;n<=num;n++)
      nnzs[findex(num,z,n)]=0;
  
 
  //Create list of sides in zone
  for (s=1;s<=kksl;s++)
    if(kstyp[findex1(s)]>0)  {   
      z=kksz[findex1(s)];
 
      for (i=1;i<=numz[findex1(z)];i++)
	if (nnzs[findex(num,z,i)]==s)
	  break;
	else if (nnzs[findex(num,z,i)]==0) {
	  nnzs[findex(num,z,i)]=s;
	  break;
	}
    } 

  kkss5= new int[kksll];
  for (s=1;s<=kksll;s++)
    kkss5[findex1(s)]=0;

  //Match adjacent side in zone
  for (s=1;s<=kksl;s++)
    if(kstyp[findex1(s)]>0)  {  
      if (kkss5[findex1(s)]==0) {
	p1=kksp1[findex1(s)];
	p2=kksp2[findex1(s)];
	z =kksz[findex1(s)];
                  
	//Find p1 & z in the list
        for (i=1;i<=numz[findex1(z)];i++) {
	  ss5=nnzs[findex(num,z,i)];
	  pp1=kksp1[findex1(ss5)];
	  pp2=kksp2[findex1(ss5)];
	  if (pp1==p2 && pp2==p1) {
	    if (kkss5[findex1(s)]==0) {
	      kkss5[findex1(s)] =ss5;
	      kkss5[findex1(ss5)]=s;
	    }
	    break;
	  }
	}
      }
      else {
	ss5=kkss5[findex1(s)];
	if (kkss5[findex1(ss5)]!=s) printf("KKSS5a error\n");
      }

    }
      
  ksleft= new int[kksll];
  for (s=1;s<=kksll;s++)
    ksleft[findex1(s)]=0;

  //Tag left sides
  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]>0) {
      z1=kksz[findex1(s)];
      z2=kksz2_ejl[findex1(s)];
      if (kztyp[findex1(z2)]<0) 
	ksleft[findex1(s)]=1;       
      else if (z1<z2) 
	ksleft[findex1(s)]=1; 
      else
	ksleft[findex1(s)]=0;
    }
    else
      ksleft[findex1(s)]=0;

  delete[] numz;
  delete[] nnzs;

  //setup edges
  nump= new int[kkpll];
  for (p=1;p<=kkpll;p++)
    nump[findex1(p)]=0;

  //Count number of sides per point
  num=0;
  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]!=0) {  
      p1=kksp1[findex1(s)];
      p2=kksp2[findex1(s)];
      nump[findex1(p1)]=nump[findex1(p1)]+1;
      if (p1!=p2) nump[findex1(p2)]=nump[findex1(p2)]+1;
      num=MAX(num,nump[findex1(p1)]);
      num=MAX(num,nump[findex1(p2)]);
    }

  num = MAX(num,1);

  nnpp= new int[kkpll*num];
  nnpe= new int[kkpll*num];
  kkse= new int[kksll];

  for (p=1;p<=kkpll;p++)
    for (n=1;n<=num;n++) {
      nnpp[findex(num,p,n)]=0;
      nnpe[findex(num,p,n)]=0;
    }

  for (s=1;s<=kksll;s++)
    kkse[findex1(s)]=0;

  kkel=0;

  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]!=0) {

      p =kksp1[findex1(s)];
      p2=kksp2[findex1(s)];
      if (p>=p2) {
	p =p2;
	p2=kksp1[findex1(s)];
      }
      if (p!=p2) {
 
	for (i=1;i<=nump[findex1(p)];i++) {
	  if (nnpp[findex(num,p,i)]==p2) {
	    k=i;
	    break;
	  }
	  else if (nnpp[findex(num,p,i)]==0) {
	    k=i;
	    kkel=kkel+1;
	    nnpp[findex(num,p,i)]=p2;
	    nnpe[findex(num,p,i)]=kkel;
	    break;
	  }
	}
      }
      else  
	k=1;
      
      kkse[findex1(s)]=nnpe[findex(num,p,k)];
    }
    else 
      kkse[findex1(s)]=0;
 
  kkell=kkel;
  kkell = MAX(kkell, 1);

  ketyp= new int[kkell];
  for (k=1;k<=kkell;k++)
    ketyp[findex1(k)]=0;
  
  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]>0) {
      e=kkse[findex1(s)];
 
      if (ketyp[findex1(e)]>-1) {
	if (kspe[findex1(s)]==0) 
	  z=kksz2_ejl[findex1(s)];
	else         				
	  z=kksz[findex1(s)];

	ketyp[findex1(e)]=kztyp[findex1(z)];
      } 
    }

  kkep1= new int[kkell];
  kkep2= new int[kkell];

  for (e=1;e<=kkell;e++) {
    kkep1[findex1(e)]=0;
    kkep2[findex1(e)]=0;
  }

  for (s=1;s<=kksl;s++) 
    if (kstyp[findex1(s)]==1 || kstyp[findex1(s)]==-1) {
      p1=kksp1[findex1(s)];
      p2=kksp2[findex1(s)];
      e=kkse[findex1(s)];
 
      if (e!=0) {
	if (p1<p2) 
	  kkep1[findex1(e)]=p1;
	else
	  kkep1[findex1(e)]=p2;
             
      }
 
    }

  for (s=1;s<=kksl;s++) 
    if (kstyp[findex1(s)]==1 || kstyp[findex1(s)]==-1) {
      p1=kksp1[findex1(s)];
      p2=kksp2[findex1(s)];
      e=kkse[findex1(s)];
 
      if (e!=0) {
	if (p1<p2) 
	  kkep2[findex1(e)]=p2;
	else
	  kkep2[findex1(e)]=p1;
             
      }
 
    }

  delete[] nump;
  delete[] nnpp;
  delete[] nnpe;

  //setup faces
  numz=new int[kkzll];
  for (z=1;z<=kkzll;z++)
    numz[findex1(z)]=0;

  max_znumf = 0;
  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]!=0) {
      z = kksz[findex1(s)];
      numz[findex1(z)] = numz[findex1(z)]+1;
      max_znumf = MAX(max_znumf,numz[findex1(z)]);
    }

  max_znumf = MAX(max_znumf,1);

  nnzz= new int[kkzll*max_znumf];
  nnzzpe= new int[kkzll*max_znumf];
  nnzf= new int[kkzll*max_znumf];

  for (z=1;z<=kkzll;z++)
    for (n=1;n<=max_znumf;n++){
      nnzz[findex(max_znumf,z,n)]=0;
      nnzzpe[findex(max_znumf,z,n)]=0;
      nnzf[findex(max_znumf,z,n)]=0;
    }

  kstag= new int[kksll];
  kksf= new int[kksll];

  f=0;

  for (s=1;s<=kksl;s++)
    if(kstyp[findex1(s)]!=0 && kstag[findex1(s)]==0) {
      f = f + 1;
      s_now = s;
      s_next = 0;
      while (s_next!=s) {
	s_next = kkss4[findex1(s_now)];
	s2 = kkss2[findex1(s_now)];
	kstag[findex1(s_now)] = 1;
	kstag[findex1(s2)]    = 1;
	kksf[findex1(s_now)] = f;
	if(kstyp[findex1(s2)]!=0) kksf[findex1(s2)] = f;
	s_now = s_next;
      }
    }

  kkfl  = f;
  kkfll = f;
      
  kftyp=new int[kkfll];
  for (f=1;f<=kkfll;f++)
    kftyp[findex1(f)]=0;

  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]>0){
      if (kspe[findex1(s)]==0)
	z=kksz2_ejl[findex1(s)];
      else
	z=kksz[findex1(s)];
      f=kksf[findex1(s)];
      if (f>0 && f<=kkfl) kftyp[findex1(f)]=kztyp[findex1(z)];
    }
  
  kfpe= new int[kkfll];
  for (f=1;f<=kkfll;f++)
    kfpe[findex1(f)]=0;  

  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]>0)
      if (kspe[findex1(s)]==1) {
	f=kksf[findex1(s)];
	kfpe[findex1(f)]=1;
      }

  kkfz2_ejl=new int[kkfll];
  kkfz2_ejlpe=new int[kkfll];
  kkfz2=new int[kkfll];
  kkfz=new int[kkfll];
  
  for (f=1;f<=kkfll;f++) {
    kkfz2_ejl[findex1(f)]=0;
    kkfz2_ejlpe[findex1(f)]=0;
    kkfz2[findex1(f)]=0;
    kkfz[findex1(f)]=0;
  }
  
  for (s=1;s<=kksl;s++)
    if (kstyp[findex1(s)]!=0){
      f=kksf[findex1(s)];
      z1=kksz[findex1(s)];
      z2=kksz2_ejl[findex1(s)];
      if (kstyp[findex1(s)]<0) z2=kksz[findex1(s)];
    
      if (f>0 && f<=kkfl) {
	if (kztyp[findex1(z2)]<0) 
	  kkfz2_ejl[findex1(f)]=z2;    
	else if (z1<z2) 
	  kkfz2_ejl[findex1(f)]=z2;
      }
    }
  
  for (s=1;s<=kksl;s++){
    z1=kksz[findex1(s)];
    z2=kksz2_ejl[findex1(s)];
    f=kksf[findex1(s)];

    if (f>0 && f<=kkfl) {
      if (kztyp[findex1(z2)]<0) 
	kkfz[findex1(f)]=z1;  
      else if (z1<z2)
	kkfz[findex1(f)]=z1;   
      else if (z2>0) 
	kkfz[findex1(f)]=z2;
    }
  }
  
  for (f=1;f<=kkfll;f++)
    kkfz2[findex1(f)]=-1;

  for (s=1;s<=kksl;s++) {
    f =kksf[findex1(s)];
    if (f>0 && f<=kkfl) {
      z1 = kkfz[findex1(f)];
      z2 = kksz2[findex1(s)];
      if (kstyp[findex1(s)]<0) z2 = kksz[findex1(s)];
      if(z1!=z2) kkfz2[findex1(f)] = z2;
    }
  }

  for (f=1;f<=kkfll;f++)    
    if(kftyp[findex1(f)]!=0 && kkfz2[findex1(f)]==-1)  printf("Error setting KKFZ2\n");

  delete[] numz;
  delete[] kstag;
  delete[] nnzz;
  delete[] nnzzpe;
  delete[] nnzf;

  //setup corners
  numz=new int[kkzll];
  for (z=1;z<=kkzll;z++)
    numz[findex1(z)]=0;

  s=kkzl;
  num=0;
  for (s=1;s<=kksl;s++) 
    if (kstyp[findex1(s)]!=0){
      z=kksz[findex1(s)];
      numz[findex1(z)]=numz[findex1(z)]+2;
      num=MAX(num,numz[findex1(z)]);
    }
  
  num = MAX(num, 1);

  nnzp= new int[num*kkzll];
  nnzc= new int[num*kkzll];

  for (z=1;z<=kkzll;z++)
    for (n=1;n<=num;n++) {
      nnzp[findex(num,z,n)]=0;
      nnzc[findex(num,z,n)]=0;
    }

  
  kksc1 = new int[kksll];
  kksc2 = new int[kksll];

  for (s=1;s<=kksll;s++) {
    kksc1[findex1(s)]=0;
    kksc2[findex1(s)]=0;
  }

  c=0;
  for (s=1;s<=kksl;s++) 
    if (kstyp[findex1(s)]!=0) {
      z=kksz[findex1(s)];
      p1=kksp1[findex1(s)];
      p2=kksp2[findex1(s)];
 
      for (i=1;i<=numz[findex1(z)];i++) 
	if (nnzp[findex(num,z,i)]==p1)  {
	  k=i;
	  break;
	}
	else if (nnzp[findex(num,z,i)]==0) {
	  k=i;
	  c=c+1;
	  nnzc[findex(num,z,i)]=c;
	  nnzp[findex(num,z,i)]=p1;
	  break;
	}

      kksc1[findex1(s)]=nnzc[findex(num,z,k)];
      for (i=1;i<=numz[findex1(z)];i++) 
	if (nnzp[findex(num,z,i)]==p2){
	  k=i;
	  break;
	}
	else if (nnzp[findex(num,z,i)]==0) {
	  k=i;
	  c=c+1;
	  nnzc[findex(num,z,i)]=c;
	  nnzp[findex(num,z,i)]=p2;
	  break;
	}
      kksc2[findex1(s)]=nnzc[findex(num,z,k)];

    }
 
 
  kkcl=c;
  kkcll=MAX(1,kkcl);
  
  kkcz=new int[kkcll];
  kkcp=new int[kkcll];
  kctyp=new int[kkcll];
  kkcs1=new int[kkcll];
  kkcs2=new int[kkcll];

  for (c=1;c<=kkcll;c++) {
    kkcz[findex1(c)]=0;
    kkcp[findex1(c)]=0;
    kctyp[findex1(c)]=0;
    kkcs1[findex1(c)]=0;
    kkcs2[findex1(c)]=0;
  }

  for (s=1;s<=kksl;s++) {
 
    c1 = kksc1[findex1(s)];
    c2 = kksc2[findex1(s)];
 
    kkcs1[findex1(c1)] = s;
    kkcs2[findex1(c2)] = s;
  }
 
  
  for (z=1;z<=kkzl;z++) 
    if (kztyp[findex1(z)]!=0) 
      for (i=1;i<=numz[findex1(z)];i++) {
	c=nnzc[findex(num,z,i)];
	if (c!=0) {
	  kkcz[findex1(c)]=z;
	  kkcp[findex1(c)]=nnzp[findex(num,z,i)];
	  kctyp[findex1(c)]=kztyp[findex1(z)];
	}
      }

  delete[] numz;
  delete[] nnzp;
  delete[] nnzc;

  num=0;
  for (p=1;p<=kkpl;p++)
    if(kptyp[findex1(p)]!=0) num++;
  printf("number of points = %d\n",num);

  num=0;
  for (p=1;p<=kkzl;p++)
    if(kztyp[findex1(p)]!=0) num++;
  printf("number of zones = %d\n",num);

  num=0;
  for (p=1;p<=kksl;p++)
    if(kstyp[findex1(p)]!=0) num++;
  printf("number of sides = %d\n",num);

  num=0;
  for (p=1;p<=kkel;p++)
    if(ketyp[findex1(p)]!=0) num++;
  printf("number of edges = %d\n",num);

  num=0;
  for (p=1;p<=kkfl;p++)
    if(kftyp[findex1(p)]!=0) num++;
  printf("number of faces = %d\n",num);  
  

}
