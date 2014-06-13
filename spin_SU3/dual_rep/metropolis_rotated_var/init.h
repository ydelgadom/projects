#ifndef _UTILS_H
#define _UTILS_H

/************************************************************************

	Subroutines in this file:

	void read_params()
	void init_fields()
	void mk_arrays()
	void rm_arrays()
	void init_vtau()
	void init_vemu()
	void init_vkappa()

**************************************************************************/

#ifdef MU
void read_params()
{
  char dummy[30];

  file.open("metro_su3_mu.start", ios::in);
  file >> dummy >> par0;
  printf(" mu0   = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dmu   = %f\n", dpar);
  file >> dummy >> npar;
  printf(" nmu   = %d\n", npar);
  file >> dummy >> kappa;
  printf(" kappa = %f\n", kappa);
  file >> dummy >> tau;
  printf(" tau   = %f\n", tau);
  file >> dummy >> nequi;
  printf(" nequi = %d\n", nequi);
  file >> dummy >> nmeas;
  printf(" nmeas = %d\n", nmeas);
  file >> dummy >> nskip;
  printf(" nskip = %d\n", nskip);
  file >> dummy >> iseed;
  printf(" iseed = %d\n", iseed);
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

//_________________________________________________________________________
#else
void read_params()
{
  char dummy[30];

  file.open("metro_su3_tau.start", ios::in);
  file >> dummy >> par0;
  printf(" tau0  = %f\n", par0);
  file >> dummy >> dpar;
  printf(" dtau  = %f\n", dpar);
  file >> dummy >> npar;
  printf(" ntau  = %d\n", npar);
  file >> dummy >> kappa;
  printf(" kappa = %f\n", kappa);
  file >> dummy >> mu;
  printf(" mu    = %f\n", mu);
  file >> dummy >> nequi;
  printf(" nequi = %d\n", nequi);
  file >> dummy >> nmeas;
  printf(" nmeas = %d\n", nmeas);
  file >> dummy >> nskip;
  printf(" nskip = %d\n", nskip);
  file >> dummy >> iseed;
  printf(" iseed = %d\n", iseed);
  file >> dummy >> outfile;
  printf(" outfile = %s\n", outfile);
  file.close();

#ifdef KAPPA0
	if ( kappa>0.0 )
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
  init_lattice();

  ndim    = 0;
  nmon[0] = 0;
  nmon[1] = 0;
  nmon[2] = 0;

  for( int i=0; i<nsite; i++ ){
    for ( int j=0; j<2; j++){
      mon[i][j] = 0;
      dim[i][2*j]   = 0;
      dim[i][2*j+1] = 0;
      dim[i][2*j+2] = 0;
      mnx[i][j] = 0;
    }
  }

  //initialize Tmn
  for (int i=0; i<MAXFLUX; i++){
    facn[i] = fact(i);
    for (int j=i; j<MAXFLUX; j++){
      Tmn[i][j] = Inm(i,j);
      Tmn[j][i] = Tmn[i][j];
    }
  }
}

//_______________________________________________________________________
void mk_arrays()
{
  nblock  = new int[4*nmeas]; 
}

//_______________________________________________________________________
  void rm_arrays()
{
  delete[] nblock;
}

//_______________________________________________________________________
void init_vtau()
{
  vtau[7]  = 1.0;
  vtau[6]  = 1.0/tau;
  vtau[5]  = vtau[6]/tau;
  vtau[4]  = vtau[5]/tau;
  vtau[3]  = vtau[4]/tau;
  vtau[2]  = vtau[3]/tau;
  vtau[1]  = vtau[2]/tau;
  vtau[0]  = vtau[1]/tau;
  vtau[8]  = tau;
  vtau[9]  = vtau[8]*tau;
  vtau[10] = vtau[9]*tau;
  vtau[11] = vtau[10]*tau;
  vtau[12] = vtau[11]*tau;
  vtau[13] = vtau[12]*tau;
  vtau[14] = vtau[13]*tau;
}

//_______________________________________________________________________
void init_vemu()
{
  emu = exp(mu);

  vemu[7]  = 1.0;
  vemu[6]  = 1.0/emu;
  vemu[5]  = vemu[6]/emu;
  vemu[4]  = vemu[5]/emu;
  vemu[3]  = vemu[4]/emu;
  vemu[2]  = vemu[3]/emu;
  vemu[1]  = vemu[2]/emu;
  vemu[0]  = vemu[1]/emu;
  vemu[8]  = emu;
  vemu[9]  = vemu[8]*emu;
  vemu[10] = vemu[9]*emu;
  vemu[11] = vemu[10]*emu;
  vemu[12] = vemu[11]*emu;
  vemu[13] = vemu[12]*emu;
  vemu[14] = vemu[13]*emu;
}

//_______________________________________________________________________
void init_vkappa()
{
  vkappa[7]  = 1.0;
  vkappa[6]  = 1.0/kappa;
  vkappa[5]  = vkappa[6]/kappa;
  vkappa[4]  = vkappa[5]/kappa;
  vkappa[3]  = vkappa[4]/kappa;
  vkappa[2]  = vkappa[3]/kappa;
  vkappa[1]  = vkappa[2]/kappa;
  vkappa[0]  = vkappa[1]/kappa;
  vkappa[8]  = kappa;
  vkappa[9]  = vkappa[8]*kappa;
  vkappa[10] = vkappa[9]*kappa;
  vkappa[11] = vkappa[10]*kappa;
  vkappa[12] = vkappa[11]*kappa;
  vkappa[13] = vkappa[12]*kappa;
  vkappa[14] = vkappa[13]*kappa;
}

#endif
