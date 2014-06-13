/***********************************************************************************
* Analysis program:
* Calculates bulk observables and susceptibilities of the Z_3 Gauge-Higgs model.
*
*	To execute: ./bin/anal_$PAR.x -f file_with_configs
*
* Compile with: make PAR=BETA anal
*
*************************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <ctime>
#include <getopt.h>
using namespace std;

fstream file;
char infile[150],outfile[150];

// Max. dimension of arrays
const int nparmax = 50;
const int nmeasmax = 1000000;

// Lattice size
int leng,leng_t,nsite;

// Global MC parameters
int nbeta,nmeas;
double beta,mu,gama;
double par[nparmax];
double dpar, par0;
int npar;

// Observables per measurement
double nlink[nmeasmax],nplaq[nmeasmax],nlink4[nmeasmax];

// Linear contribution to the second derivatives
double ulin[nmeasmax],nlin[nmeasmax];

// Weights
double mm[3],dmm[3],d2mm[3],rr[3],qq[3];
double bb,bbt,dbbt,aa,kk;

// Average observables
double linkaver[nparmax],linkerr[nparmax]; // spatial link occupation number
double plaqaver[nparmax],plaqerr[nparmax]; // plaquette occup. number
double linkaver4[nparmax],linkerr4[nparmax]; // temporal link occup. number
double slinkaver[nparmax],slinkerr[nparmax]; // susceptibility of spatial link
double slinkaver4[nparmax],slinkerr4[nparmax]; //  susceptibility of temporal link
double splaqaver[nparmax],splaqerr[nparmax]; // susceptibility of plaquette


// Subroutines
void processdata( int &argc, char *argv[] );
void finalize(int it);
void calculate_bbweights();
void calculate_bbtweights();
void calculate_link4weights();

char texthelp[]="Usage: exec -f [FILE]\n"
		"Analysis program of the Z(3) gauge-Higgs model in the dual rep.\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -f, --F  File with the configurations\n"
		"Report bugs to ydelgado83@gmail.com\n";

//-----------------------------------------------------------------------------   
int main( int argc, char *argv[] )
{
  printf("\nProgram analysis.cpp\n");

  processdata( argc, argv );
   
	// print observables in output file
  file.open(outfile, ios::out | ios::trunc );
  for (int ipar = 0; ipar<=npar ; ipar++ )
	{  
    file << par[ipar]        << " "
         << linkaver[ipar]   << " " << linkerr[ipar]   << " "
         << slinkaver[ipar]  << " " << slinkerr[ipar]  << " "
         << linkaver4[ipar]  << " " << linkerr4[ipar]  << " "
         << slinkaver4[ipar] << " " << slinkerr4[ipar] << " "
         << plaqaver[ipar]   << " " << plaqerr[ipar]   << " "
         << splaqaver[ipar]  << " " << splaqerr[ipar]  << " "
         << endl;
  }
  file.close();
  printf("\nDone.\n");

}
  
//!------------------------------------------------------------------------
//!------------------------------------------------------------------------
void finalize(int it)
{
	/*
		This subrouting computes the average values of the
		observables and the errors (Jackknife is used
		for the susceptibilities)
	*/
  int    im;
  double esum,efsum,efblock,cfsum;
  double msum,mfsum,mfblock,sfsum;
  double esum4,efsum4,efblock4,cfsum4;

  msum = 0.0;
  esum = 0.0;
  esum4 = 0.0;
  for(im=0 ; im<nmeas ; im++){
    esum = esum + nlink[im];
    esum4 = esum4 + nlink4[im];
    msum = msum + nplaq[im];
  }
  linkaver[it] = esum/(double)nmeas;
  linkaver4[it] = esum4/(double)nmeas;
  plaqaver[it] = msum/(double)nmeas;

  efsum = 0.0;
  efsum4 = 0.0;
  mfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    efsum += ( linkaver[it] - nlink[im] )*( linkaver[it] - nlink[im] );
    efsum4 += ( linkaver4[it] - nlink4[im] )*( linkaver4[it] - nlink4[im] );
    mfsum += ( plaqaver[it] - nplaq[im] )*( plaqaver[it] - nplaq[im] );
  }
  linkerr[it] = sqrt(efsum)/(double)nmeas;
  linkerr4[it] = sqrt(efsum4)/(double)nmeas;
  plaqerr[it] = sqrt(mfsum)/(double)nmeas;

  for(im=0 ; im<nmeas ; im++){
    efsum  += 0.;       // no linear term for susceptibility
    efsum4 += nlin[im]; // n
    mfsum  += ulin[im]; // energy
  }
  slinkaver[it] = efsum/(double)nmeas;
  slinkaver4[it] = efsum4/(double)nmeas;
  splaqaver[it] = mfsum/(double)nmeas;

  cfsum  = 0.0;
  cfsum4  = 0.0;
  sfsum  = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock  = ( efsum - ( linkaver[it] - nlink[im] )*( linkaver[it] - nlink[im] ) - 0. )/double(nmeas-1);
    efblock4 = ( efsum4 - ( linkaver4[it] - nlink4[im] )*( linkaver4[it] - nlink4[im] ) - nlin[im] )/double(nmeas-1);
    mfblock  = ( mfsum - ( plaqaver[it] - nplaq[im] )*( plaqaver[it] - nplaq[im] ) - ulin[im] )/double(nmeas-1);

    cfsum  += ( slinkaver[it] - efblock )*( slinkaver[it] - efblock );
    cfsum4 += ( slinkaver4[it] - efblock4 )*( slinkaver4[it] - efblock4 );
    sfsum  += ( splaqaver[it] - mfblock )*( splaqaver[it] - mfblock );
  }  
  slinkerr[it]  = sqrt(cfsum);
  slinkerr4[it] = sqrt(cfsum4);
  splaqerr[it]  = sqrt(sfsum);
   
  linkaver[it]  /= (double)nsite;
  linkerr[it]   /= (double)nsite;
  slinkaver[it] /= (double)nsite;
  slinkerr[it]  /= (double)nsite;

  linkaver4[it]  /= (double)nsite;
  linkerr4[it]   /= (double)nsite;
  slinkaver4[it] /= (double)nsite;
  slinkerr4[it]  /= (double)nsite;

  plaqaver[it]  /= (double)nsite;
  plaqaver[it]  = plaqaver[it]/6. + bbt;
  plaqerr[it]   /= (double)nsite;
  plaqerr[it]   /= 6.;
  splaqaver[it] /= (double)nsite;
  splaqaver[it] += 6.*dbbt;
  splaqaver[it] /= 6.;
  splaqerr[it]  /= (double)nsite;
  splaqerr[it]  /= 6.;
}

//!------------------------------------------------------------------------
inline void getfilename( int &argc, char *argv[] )
{
	// Read name of file with configurations
	if(argc<1) cout << endl << texthelp << endl;

  int c;    
	while (1)
	{
		static struct option long_options[] =
			{
			/* These options don't set a flag.
			We distinguish them by their indices. */
			{"F", required_argument, 0, 'f'},
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "f:",
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c){
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");

			case 'f':
				sprintf(infile,"%s",optarg);
			  sprintf(outfile,"%sobs",infile);
				cout << "Input File: " << infile << endl;
  			cout << "Output file: " << outfile << endl;
				break;

			default:
				cout << endl << texthelp << endl;
		}
	} 
}
   
//!------------------------------------------------------------------------
void processdata( int &argc, char *argv[] )
{
	// Local MC parameters
  int nequi,nskip,iseed;

	// Occupation numbers
  int nnlink[3],nnplaq[3],nnlink4[3];

	// Array with measurements from gen.cpp
  int *nblock;
  
	getfilename( argc, argv );

  file.open(infile, ios::in | ios::binary);
  file.read((char*)&leng,sizeof(int));
  file.read((char*)&leng_t,sizeof(int));
  file.read((char*)&gama,sizeof(double));
#ifdef MU
  file.read((char*)&beta,sizeof(double));
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
#else
  file.read((char*)&mu,sizeof(double));
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
#endif
  file.read((char*)&nequi,sizeof(int));
  file.read((char*)&nmeas,sizeof(int));
  file.read((char*)&nskip,sizeof(int));
  file.read((char*)&iseed,sizeof(int));
     
  printf(" leng   = %d\n", leng);
  printf(" leng_t = %d\n", leng_t);
  printf(" gamma  = %f\n", gama);
#ifdef MU
  printf(" beta   = %f\n", beta);
  printf(" mu0    = %f\n", par0);
  printf(" dmu    = %f\n", dpar);
  printf(" nmu    = %d\n", npar);
#else
  printf(" mu     = %f\n", mu);
  printf(" beta0  = %f\n", par0);
  printf(" dbeta  = %f\n", dpar);
  printf(" nbeta  = %d\n", npar);
#endif
  printf(" nequi  = %d\n", nequi);
  printf(" nmeas  = %d\n", nmeas);
  printf(" nskip  = %d\n", nskip);
  printf(" iseed  = %d\n", iseed);
  
  nsite  = leng*leng*leng*leng_t;
  nblock = new int[9*nmeas];

  calculate_bbweights();

  for( int ipar=0 ; ipar<=npar ; ipar++ )
  {    
    par[ipar] = par0 + ipar*dpar;
#ifdef MU
		mu = par[ipar];
#else
		beta = par[ipar];
#endif
    cout << "Reading beta= " << beta << " - mu = " << mu << endl;

		// Compute weights
    calculate_bbtweights();
  	calculate_link4weights();
	  		   
    file.read( (char*)nblock, 9*nmeas*sizeof(int) );

    for ( int imeas = 0; imeas<nmeas; imeas++ )
    {
       nnlink[0] = nblock[imeas*9+0];
       nnlink[1] = nblock[imeas*9+1];
       nnlink[2] = nblock[imeas*9+2];

       nnlink4[0] = nblock[imeas*9+3];
       nnlink4[1] = nblock[imeas*9+4];
       nnlink4[2] = nblock[imeas*9+5];

       nnplaq[0] = nblock[imeas*9+6];
       nnplaq[1] = nblock[imeas*9+7];
       nnplaq[2] = nblock[imeas*9+8];
   
       nlink[imeas]  = nnlink[0] + nnlink[2];

       //nlink4[imeas] = nnlink4[0] + nnlink4[2];
       //nlin[imeas]  = 0.;
       nlink4[imeas] = rr[0]*nnlink4[0] + rr[1]*nnlink[1] + rr[2]*nnlink4[2];  // n
       nlin[imeas]   = qq[0]*nnlink4[0] + qq[1]*nnlink[1] + qq[2]*nnlink4[2];

       //nplaq[imeas]  = nnplaq[0] + nnplaq[2];
       //ulin[imeas] = 0.;
       nplaq[imeas] = aa*(nnplaq[0] + nnplaq[2]); // U
       ulin[imeas]  = kk*(nnplaq[0] + nnplaq[2]);
    }
    finalize(ipar);
    
  }
  
  file.close();
  //delete[] nblock;
   
}

//_________________________________________________________________________
void calculate_bbweights()
{
	bb  = ( exp(2.*gama) - exp(-gama) )/( exp(2.*gama) + 2.*exp(-gama) );
}

//_________________________________________________________________________
void calculate_bbtweights()
{
  bbt = ( exp(beta) - exp(-beta/2.) )/( exp(beta) + 2.*exp(-beta/2.) );
  dbbt = 9. / ( 2.*( exp(1.5*beta) + 4.*exp(-1.5*beta) + 4. ) ) ; 
  aa  = 9. / ( 2.*( exp(1.5*beta) - 2.*exp(-1.5*beta) + 1. ) ) ;  // dbbt/bbt
  kk  =  -4.5 * ( (1.5*exp(1.5*beta) - 6.*exp(-1.5*beta))
			/pow((exp(1.5*beta) + 4.*exp(-1.5*beta) + 4.),2) ) / bbt - aa*aa;
}

//_________________________________________________________________________
void calculate_link4weights()
{
  int i;
  const double pi2over3 = 2.0*3.141592653589793/3.0;
  double sinhu = sinh(mu);
  double coshu = cosh(mu);
  double sqrt3 = sqrt(3.);

  for (i = -1; i<=1; i++)
  {
    mm[i+1] = ( exp(2.*gama*coshu) + 
                2.*exp(-gama*coshu)*cos(sinhu*sqrt3*gama - i*pi2over3) ) / 3.0;

    dmm[i+1] = 2.* gama * ( sinhu*exp(2.*gama*coshu) - 
                            exp(-gama*coshu) * sinhu * cos(sinhu*sqrt3*gama - i*pi2over3) -
                            exp(-gama*coshu) * sqrt3 * coshu * sin(sinhu*sqrt3*gama - i*pi2over3)  ) / 3.0;

    d2mm[i+1] = 2.* gama * ( exp(2.*gama*coshu) * (2.*gama*sinhu*sinhu + coshu) +
                             exp(-gama*coshu) * (gama*sinhu*sinhu - coshu - 3.*gama*coshu*coshu) * cos(sinhu*sqrt3*gama - i*pi2over3) +
                             exp(-gama*coshu) * sqrt3 * sinhu * (2.*gama*coshu - 1.) * sin(sinhu*sqrt3*gama - i*pi2over3)  ) / 3.0;

    rr[i+1] = dmm[i+1]/mm[i+1];

    qq[i+1] = d2mm[i+1]/mm[i+1] - rr[i+1]*rr[i+1];

    cout << "M["<< i <<"] = " << mm[i+1] << endl;
  }
  cout << endl;

}
