#ifndef _INIT_H
#define _INIT_H

/************************************************************************

  Subroutines in this file:

  void read_params()
  void init_fields()
  void mk_arrays()
  void rm_arrays()
  void init_vtau()
  void init_veta()
  void init_vetabar()

**************************************************************************/

#ifdef MU
void read_params()
{
  char dummy[30];

  file.open("metro_su3_mu.start", ios::in);
  
  file >> dummy >> par0;
  printf(" mu0    = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dmu    = %f\n", dpar);
  file >> dummy >> npar;
  printf(" nmu    = %d\n", npar);
  file >> dummy >> kappa;
  printf(" kappa  = %f\n", kappa);
  file >> dummy >> tau;
  printf(" tau    = %f\n", tau);
  file >> dummy >> nequi;
  printf(" nequi  = %d\n", nequi);
  file >> dummy >> nmeas;
  printf(" nmeas  = %d\n", nmeas);
  file >> dummy >> nskip;
  printf(" nskip  = %d\n", nskip);
  file >> dummy >> iseed;
  printf(" iseed  = %d\n", iseed);
  file >> dummy >> outfile;
  printf(" outfile = %s\n", outfile);
  
  file.close();
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&par0,sizeof(double));
  file.write((char*)&dpar,sizeof(double));
  file.write((char*)&npar,sizeof(int));
  file.write((char*)&kappa,sizeof(double));
  file.write((char*)&tau,sizeof(double));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));
}
#endif

//_________________________________________________________________________
#ifndef MU
void read_params()
{
  char dummy[30];

  file.open("metro_su3_tau.start", ios::in);
  
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

#ifdef KAPPA0
  if (kappa>0)
  {
    cout << "ERROR: kappa !neq 0" << endl;
    cout << "Setting kappa to 0" << endl;
    kappa = 0.0;  
  }
#endif
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
 
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

  ndim    = 0;
  nmon[0] = 0;
  nmon[1] = 0;

  for( i=0; i<nsite; i++ ){
    for (j=0; j<2; j++){
      mon[i][j] = 0;
      dim[i][2*j]   = 0;
      dim[i][2*j+1] = 0;
      dim[i][2*j+2] = 0;
      mnx[i][j] = 0;
    }
  }

  //initialize Tmn
  for (i=0; i<MAXFLUX; i++){
    facn[i] = fact(i);
    for (j=i; j<MAXFLUX; j++){
      Tmn[i][j] = Inm(i,j);
      Tmn[j][i] = Tmn[i][j];
    }
  }
}

//_______________________________________________________________________
void mk_arrays()
{
  init_lattice();

  nblock  = new int[3*nmeas]; 
}


//_______________________________________________________________________
void rm_arrays()
{
  delete[] nblock;
}

//_______________________________________________________________________
void init_vtau()
{
  vtau[3] = 1.0;
  vtau[2] = 1.0/tau;
  vtau[1] = vtau[2]/tau;
  vtau[0] = vtau[1]/tau;
  vtau[4] = tau;
  vtau[5] = vtau[4]*tau;
  vtau[6] = vtau[5]*tau;
}

//_______________________________________________________________________
void init_veta()
{
  veta[3] = 1.0;
  veta[2] = 1.0/eta;
  veta[1] = veta[2]/eta;
  veta[0] = veta[1]/eta;
  veta[4] = eta;
  veta[5] = veta[4]*eta;
  veta[6] = veta[5]*eta;
}

//_______________________________________________________________________
void init_vetabar()
{
  vetabar[3] = 1.0;
  vetabar[2] = 1.0/etabar;
  vetabar[1] = vetabar[2]/etabar;
  vetabar[0] = vetabar[1]/etabar;
  vetabar[4] = etabar;
  vetabar[5] = vetabar[4]*etabar;
  vetabar[6] = vetabar[5]*etabar;
}

#endif
