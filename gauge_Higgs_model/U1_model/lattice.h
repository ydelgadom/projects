#ifndef _LATTICE_H
#define _LATTICE_H

// Variables
const int leng   = SIZE;
const int leng_t = leng;
const int nsite = leng*leng*leng*leng_t;

int neib[nsite][8];


// Subroutines
void init_lattice()
{
  int leng2 = leng*leng;
  int leng3 = leng2*leng;

  int i1,i2,i3,i4;
  for (i1 = 0; i1<leng ; i1++ ){
    int i1p = i1 + 1;
    int i1m = i1 - 1;
    if (i1 == (leng-1) ) i1p = 0;
    if (i1 == 0)         i1m = leng-1;

  for (i2 = 0; i2<leng ; i2++){
    int i2p = i2 + 1;
    int i2m = i2 - 1;
    if (i2 == (leng-1) ) i2p = 0;
    if (i2 == 0)         i2m = leng-1;

  for (i3 = 0; i3<leng ; i3++){
    int i3p = i3 + 1;
    int i3m = i3 - 1;
    if (i3 == (leng-1) ) i3p = 0;
    if (i3 == 0)         i3m = leng-1;

  for (i4 = 0; i4<leng_t ; i4++){
    int i4p = i4 + 1;
    int i4m = i4 - 1;
    if (i4 == (leng_t-1) ) i4p = 0;
    if (i4 == 0)           i4m = leng_t-1;

  int is = i1 + i2*leng + i3*leng2 + i4*leng3;  

  int isp1 = i1p + i2*leng  + i3*leng2  + i4*leng3;   
  int isp2 = i1  + i2p*leng + i3*leng2  + i4*leng3; 
  int isp3 = i1  + i2*leng  + i3p*leng2 + i4*leng3;
  int isp4 = i1  + i2*leng  + i3*leng2  + i4p*leng3;
 
  int ism1 = i1m + i2*leng  + i3*leng2  + i4*leng3;
  int ism2 = i1  + i2m*leng + i3*leng2  + i4*leng3;
  int ism3 = i1  + i2*leng  + i3m*leng2 + i4*leng3;
  int ism4 = i1  + i2*leng  + i3*leng2  + i4m*leng3;

  // fill in the neighborhood array

  neib[is][0] = isp1;
  neib[is][1] = isp2;
  neib[is][2] = isp3;
  neib[is][3] = isp4;
 
  neib[is][4] = ism1;
  neib[is][5] = ism2;
  neib[is][6] = ism3;
  neib[is][7] = ism4;

  }
  }
  }
  }

}

#endif
