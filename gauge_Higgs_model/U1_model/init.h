#ifndef _INIT_H
#define _INIT_H

void read_params( )
{
  char dummy[30];

  file.open("worm_beta.start", ios::in);
  
  file >> dummy >> kappa;
  file >> dummy >> lambda;
  file >> dummy >> beta0;
  file >> dummy >> dbeta;
  file >> dummy >> nbeta;
  file >> dummy >> nequi;
  file >> dummy >> nmeas;
  file >> dummy >> nskip;
  file >> dummy >> iseed;
  file >> dummy >> outfile;
  
  file.close();

  printf(" kappa   = %f\n", kappa);
  printf(" lambda  = %f\n", lambda);
  printf(" beta0   = %f\n", beta0);
  printf(" dbeta   = %f\n", dbeta);
  printf(" nbeta   = %d\n", nbeta);
  printf(" nequi   = %d\n", nequi);
  printf(" nmeas   = %d\n", nmeas);
  printf(" nskip   = %d\n", nskip);
  printf(" iseed   = %d\n", iseed);
  printf(" outfile = %s\n", outfile);
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&leng_t,sizeof(int));
  file.write((char*)&kappa,sizeof(double));
  file.write((char*)&lambda,sizeof(double));
  file.write((char*)&beta0,sizeof(double));
  file.write((char*)&dbeta,sizeof(double));
  file.write((char*)&nbeta,sizeof(int));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));
}

//_______________________________________________________________________
void init( )
{
  read_params( );
  rlxd_init( 1, iseed );
  init_lattice( );
  
  memset ( nblock, 0, (2*LENGTH+LENGTH_FLUX)*sizeof(int) );
  memset ( nplaq,  0, LENGTH*sizeof(int) );
  memset ( nflux,  0, LENGTH_FLUX*sizeof(int) );
  memset ( nlink,  0, LENGTH*sizeof(int) );

  nplaq[0] = 6*nsite;
  nflux[0] = nsite;
  nlink[0] = 4*nsite;

  memset ( vlink,  0, nsite*4*sizeof(int) );
  memset ( vllink, 0, nsite*4*sizeof(int) );
  memset ( vplaq,  0, nsite*6*sizeof(int) );
  memset ( vflux,  0, nsite*sizeof(int) );

  int i;
  for ( i=0; i<LENGTH; i++ )
    fac[i] = log( fact( i ) );
}

#endif
