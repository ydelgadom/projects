#ifndef _MULTIGRID_H
#define _MULTIGRID_H

  /*
    Performs 90% of equilibration steps 
    on a smaller part of the lattice: sub_leng^3*Nt (leng must be a multiple of sub_leng).
    Then, this sublattice is replicated in the full lattice.
  */

int multigrid_equi( )
{
  // Array with neighbours in the smaller lattice
  int sneib[snsite][8];
  init_lattice( sub_leng, sneib );

  // Do equilibration sweeps
  if ( nsweeps( 9*nequi/10, sub_leng, snsite, sneib ) < 0 ) return -1;
#ifdef CHECK
  check( snsite, sneib );
#endif
  int tmpplaq[snsite][6], tmpflux[snsite][2], tmplink[snsite][8], tmpllink[snsite][8];

  memcpy( &tmpplaq[0][0], &vplaq[0][0], snsite*6*sizeof(int) );
  memcpy( &tmpflux[0][0], &vflux[0][0], snsite*2*sizeof(int) );
  memcpy( &tmplink[0][0], &vlink[0][0], snsite*8*sizeof(int) );
  memcpy( &tmpllink[0][0], &vllink[0][0], snsite*8*sizeof(int) );

  // Replicate sublattice
  int factor = (leng/sub_leng);
  int il[3];
  for ( il[2]=0 ; il[2]<factor ; il[2]++ )
  for ( il[1]=0 ; il[1]<factor ; il[1]++ )
  for ( il[0]=0 ; il[0]<factor ; il[0]++ )
  {
    int factor0 = ( il[0] + (il[1] + il[2]*leng)*leng ) * sub_leng ;
    int i[4];
    for ( i[3]=0 ; i[3]<leng_t ; i[3]++ )
    for ( i[2]=0 ; i[2]<sub_leng ; i[2]++ )
    for ( i[1]=0 ; i[1]<sub_leng ; i[1]++ )
    for ( i[0]=0 ; i[0]<sub_leng ; i[0]++ )
    {
      int is = i[0] + sub_leng*( i[1] + sub_leng*( i[2] + i[3]*sub_leng ) );
      int is_new = i[0] + leng*( i[1] + leng*( i[2] + i[3]*leng ) ) + factor0;

      vflux[is_new][0] = tmpflux[is][0];
      vflux[is_new][1] = tmpflux[is][1];

      vplaq[is_new][0] = tmpplaq[is][0];
      vplaq[is_new][1] = tmpplaq[is][1];
      vplaq[is_new][2] = tmpplaq[is][2];
      vplaq[is_new][3] = tmpplaq[is][3];
      vplaq[is_new][4] = tmpplaq[is][4];
      vplaq[is_new][5] = tmpplaq[is][5];

      vlink[is_new][0] = tmplink[is][0];
      vlink[is_new][1] = tmplink[is][1];
      vlink[is_new][2] = tmplink[is][2];
      vlink[is_new][3] = tmplink[is][3];
      vlink[is_new][4] = tmplink[is][4];
      vlink[is_new][5] = tmplink[is][5];
      vlink[is_new][6] = tmplink[is][6];
      vlink[is_new][7] = tmplink[is][7];

      vllink[is_new][0] = tmpllink[is][0];
      vllink[is_new][1] = tmpllink[is][1];
      vllink[is_new][2] = tmpllink[is][2];
      vllink[is_new][3] = tmpllink[is][3];
      vllink[is_new][4] = tmpllink[is][4];
      vllink[is_new][5] = tmpllink[is][5];
      vllink[is_new][6] = tmpllink[is][6];
      vllink[is_new][7] = tmpllink[is][7];
    }
  }
  factor = factor*factor*factor;
  nlink[0] *= factor;
  nlink[1] *= factor;

  int i;
  for ( i=0 ; i<LENGTH ; i++ )
  {
    nflux[0][i] *= factor;
    nflux[1][i] *= factor;
    nplaq[i] *= factor;
  }
  return 1;
}

#endif
