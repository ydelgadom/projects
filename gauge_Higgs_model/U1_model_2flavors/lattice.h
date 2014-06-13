#ifndef _LATTICE_H
#define _LATTICE_H

// Lattice size
const int leng = NS;
const int leng_t = NT;
const int nsite = leng*leng*leng*leng_t;

// Dimensions of sublattice
const int sub_leng = 4;
const int snsite = sub_leng*sub_leng*sub_leng*leng_t;

// Neiborghs indices (full lattice)
int vneib[nsite][8];


//_________________________________________________________________________
void init_lattice( int vleng, int neib[][8] )
{
  int leng2 = vleng*vleng;
  int leng3 = leng2*vleng;

  int i1,i2,i3,i4;
  for (i1 = 0; i1<vleng ; ++i1 ){
    int i1p = i1 + 1;
    int i1m = i1 - 1;
    if (i1 == (vleng-1) ) i1p = 0;
    if (i1 == 0) i1m = vleng-1;

  for (i2 = 0; i2<vleng ; ++i2){
    int i2p = i2 + 1;
    int i2m = i2 - 1;
    if (i2 == (vleng-1) ) i2p = 0;
    if (i2 == 0) i2m = vleng-1;

  for (i3 = 0; i3<vleng ; ++i3){
    int i3p = i3 + 1;
    int i3m = i3 - 1;
    if (i3 == (vleng-1) ) i3p = 0;
    if (i3 == 0)  i3m = vleng-1;

  for (i4 = 0; i4<leng_t ; ++i4){
    int i4p = i4 + 1;
    int i4m = i4 - 1;
    if (i4 == (leng_t-1) ) i4p = 0;
    if (i4 == 0) i4m = leng_t-1;

  int is = i1 + i2*vleng + i3*leng2 + i4*leng3;  

  // fill in the neighborhood array

  neib[is][0] = i1p + i2*vleng  + i3*leng2  + i4*leng3;
  neib[is][1] = i1  + i2p*vleng + i3*leng2  + i4*leng3; 
  neib[is][2] = i1  + i2*vleng  + i3p*leng2 + i4*leng3;
  neib[is][3] = i1  + i2*vleng  + i3*leng2  + i4p*leng3;
 
  neib[is][4] = i1m + i2*vleng  + i3*leng2  + i4*leng3;
  neib[is][5] = i1  + i2m*vleng + i3*leng2  + i4*leng3;
  neib[is][6] = i1  + i2*vleng  + i3m*leng2 + i4*leng3;
  neib[is][7] = i1  + i2*vleng  + i3*leng2  + i4m*leng3;

  }
  }
  }
  }
}

#endif
