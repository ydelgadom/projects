#ifndef _LATTICE_H
#define _LATTICE_H

//-----------------------------
// VARIABLES
//-----------------------------

const int leng  = SIZE;
const int nsite = leng*leng*leng;
int neib[nsite][6];

//-------------------------------
// SUBROUTINES
//-------------------------------

//_________________________________________________________________________

void init_lattice()
{

  int i1,i2,i3,i1p,i2p,i3p,i1m,i2m,i3m;
  int is,isp1,isp2,isp3,ism1,ism2,ism3;

  for (i1 = 0; i1<leng ; ++i1 ){
    i1p = i1 + 1;
    i1m = i1 - 1;
    if (i1 == (leng-1) ) i1p = 0;
    if (i1 == 0)         i1m = leng-1;

  for (i2 = 0; i2<leng ; ++i2){
    i2p = i2 + 1;
    i2m = i2 - 1;
    if (i2 == (leng-1) ) i2p = 0;
    if (i2 == 0)         i2m = leng-1;

  for (i3 = 0; i3<leng ; ++i3){
    i3p = i3 + 1;
    i3m = i3 - 1;
    if (i3 == (leng-1) ) i3p = 0;
    if (i3 == 0)         i3m = leng-1;
 
  // compute the site address and the addresses of the sites shifted
  // by one unit in each direction

  is = i1 + i2*leng + i3*leng*leng;  

  isp1 = i1p + i2*leng +  i3*leng*leng;   
  isp2 = i1 +  i2p*leng + i3*leng*leng; 
  isp3 = i1 +  i2*leng +  i3p*leng*leng;  
 
  ism1 = i1m + i2*leng +  i3*leng*leng;  
  ism2 = i1 +  i2m*leng + i3*leng*leng; 
  ism3 = i1 +  i2*leng +  i3m*leng*leng;
 

  // fill in the neighborhood array

  neib[is][0] = isp1;
  neib[is][1] = isp2;
  neib[is][2] = isp3;
 
  neib[is][3] = ism1;
  neib[is][4] = ism2;
  neib[is][5] = ism3;
 
  //printf("%d - %d %d %d %d %d %d\n",is,isp1,isp2,isp3,ism1,ism2,ism3);

  }
  }
  }

}


#endif 
