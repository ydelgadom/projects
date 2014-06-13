#ifndef _INIT_H
#define _INIT_H

/***********************************************************************

	This file containes the following subroutines:

	void read_params()
	void init()
	void mk_arrays()
	void rm_arrays()

***********************************************************************/

#ifdef MU
void read_params()
{
  /* function to read input parameters from input file */

  char dummy[30];

  file.open("worm_mu.start", ios::in);
  file >> dummy >> gama;
  printf(" gamma   = %f\n", gama);  
  file >> dummy >> beta;
  printf(" beta    = %f\n", beta);
  file >> dummy >> par0;
  printf(" mu0     = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dmu     = %f\n", dpar);
  file >> dummy >> npar;
  printf(" nmu     = %d\n", npar);
  file >> dummy >> nequi;
  printf(" nequi   = %d\n", nequi);
  file >> dummy >> nmeas;
  printf(" nmeas   = %d\n", nmeas);
  file >> dummy >> nskip;
  printf(" nskip   = %d\n", nskip);
  file >> dummy >> iseed;
  printf(" iseed   = %d\n", iseed);
  file >> dummy >> outfile;
  printf(" outfile = %s\n", outfile);
  file.close();
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&leng_t,sizeof(int));
  file.write((char*)&gama,sizeof(double));
  file.write((char*)&beta,sizeof(double));
  file.write((char*)&par0,sizeof(double));
  file.write((char*)&dpar,sizeof(double));
  file.write((char*)&npar,sizeof(int));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));

}

//_________________________________________________________________________
#else
void read_params()
{
  char dummy[30];

  file.open("worm_beta.start", ios::in);
  file >> dummy >> gama;
  printf(" gamma   = %f\n", gama);  
  file >> dummy >> mu;
  printf(" mu      = %f\n", mu);
  file >> dummy >> par0;
  printf(" beta0   = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dbeta   = %f\n", dpar);
  file >> dummy >> npar;
  printf(" nbeta   = %d\n", npar);
  file >> dummy >> nequi;
  printf(" nequi   = %d\n", nequi);
  file >> dummy >> nmeas;
  printf(" nmeas   = %d\n", nmeas);
  file >> dummy >> nskip;
  printf(" nskip   = %d\n", nskip);
  file >> dummy >> iseed;
  printf(" iseed   = %d\n", iseed);
  file >> dummy >> outfile;
  printf(" outfile = %s\n", outfile);
  file.close();
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&leng_t,sizeof(int));
  file.write((char*)&gama,sizeof(double));
  file.write((char*)&mu,sizeof(double));
  file.write((char*)&par0,sizeof(double));
  file.write((char*)&dpar,sizeof(double));
  file.write((char*)&npar,sizeof(int));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));

}
#endif

//_______________________________________________________________________
void mk_arrays()
{
  nblock = new int[9*nmeas]; 
}

//_______________________________________________________________________
void rm_arrays()
{
  delete[] nblock;
}

//_________________________________________________________________________
void init()
{
  read_params();
  rlxd_init(1,iseed);
  init_lattice();
  mk_arrays();

  nlink[0][0] = 0;
  nlink[0][1] = nsite;
  nlink[0][2] = 0;
  nlink[1][0] = 0;
  nlink[1][1] = nsite;
  nlink[1][2] = 0;
  nlink[2][0] = 0;
  nlink[2][1] = nsite;
  nlink[2][2] = 0;
  nlink[3][0] = 0;
  nlink[3][1] = nsite;
  nlink[3][2] = 0;

  nplaq[0] = 0;
  nplaq[1] = 6*nsite;
  nplaq[2] = 0;

  for(int i=0; i<nsite; i++)
  {
    for(int j=0; j<4; j++)
    {
      vlink[i][j] = 1;
      vplaq[i][j] = 1;
    }
    vplaq[i][4] = 1;
    vplaq[i][5] = 1;
  }

	calculate_bb_weight( );

}

#endif
