#ifndef _INIT_H
#define _INIT_H

/**************************************************************************

	This file contains:

	void read_params( )
	void init_fields()
	void mk_arrays()
	void rm_arrays()

***************************************************************************/

//_________________________________________________________________________
#ifdef KAPPA
void read_params()
{
  char dummy[30];

  file.open("worm_kappa.start", ios::in);
  
  file >> dummy >> tau;
  printf(" tau     = %f\n", tau);
  file >> dummy >> mu;
  printf(" mu      = %f\n", mu);
  file >> dummy >> par0;
  printf(" kappa0  = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dkappa  = %f\n", dpar);
  file >> dummy >> npar;
  printf(" nkappa  = %d\n", npar);
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
   
  file.open(outfile,ios::trunc | ios::out);
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&tau,sizeof(double));
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

//_________________________________________________________________________
#ifdef MU
void read_params()
{
  char dummy[30];

  file.open("worm_mu.start", ios::in);
  
  file >> dummy >> tau;
  printf(" tau     = %f\n", tau);
  file >> dummy >> kappa;
  printf(" kappa   = %f\n", kappa);
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
   
  file.open(outfile,ios::trunc | ios::out);
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&tau,sizeof(double));
  file.write((char*)&kappa,sizeof(double));
  file.write((char*)&par0,sizeof(double));
  file.write((char*)&dpar,sizeof(double));
  file.write((char*)&npar,sizeof(int));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));

}
#endif

//_________________________________________________________________________
#ifdef TAU
void read_params()
{
  char dummy[30];

  file.open("worm_tau.start", ios::in);
  
  file >> dummy >> par0;
  printf(" tau0    = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dtau    = %f\n", dpar);
  file >> dummy >> npar;
  printf(" ntau    = %d\n", npar);
  file >> dummy >> kappa;
  printf(" kappa   = %f\n", kappa);
  file >> dummy >> mu;
  printf(" mu      = %f\n", mu);
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
 
  file.open(outfile,ios::trunc | ios::out);
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&par0,sizeof(double));
  file.write((char*)&dpar,sizeof(double));
  file.write((char*)&npar,sizeof(int));
  file.write((char*)&kappa,sizeof(double));
  file.write((char*)&mu,sizeof(double));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));
}
#endif

//_________________________________________________________________________
void init_fields()
{
  int i,j;

  init_lattice();

  nmon[0] = 0;
  nmon[1] = nsite;
  nmon[2] = 0;

  ndim[0] = 0;
  ndim[1] = 3*nsite;
  ndim[2] = 0;

  for( i=0; i<nsite; ++i ){
    mon[i] = 1; 
    for (j=0; j<3; ++j)
      dim[i][j] = 1;
  }

  triadd[0][0] = 2;
  triadd[0][1] = 0;
  triadd[0][2] = 1;
  triadd[1][0] = 0;
  triadd[1][1] = 1;
  triadd[1][2] = 2;
  triadd[2][0] = 1;  
  triadd[2][1] = 2;
  triadd[2][2] = 0;
}


//_______________________________________________________________________
void mk_arrays()
{
  nblock = new int[6*nmeas]; 
}

//_______________________________________________________________________
void rm_arrays()
{
  delete[] nblock;
}

#endif
