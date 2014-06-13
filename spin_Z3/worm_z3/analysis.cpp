/************************************************************************
* Analysis program for different values of $PAR (tau, kappa or mu)
*	in the dual representation of the Z(3) spin model.  
*	See references in ../papers.
*
*	To execute: ./bin/anal_$PAR.x -f file_with_configs
*
*	By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <ctime>
#include <getopt.h>
using namespace std;

char infile[150],outfile[150];
fstream file;

// Array dimensions
const int nparmax  = 42;
const int nmeasmax = 1000000;

// monomer and dimer weights
double monoweight[3];
double cc, cccc, bb, bbbb;

// Lattice size
int leng,nsite;

// Global MC parameters
int nmeas, npar;

// Global parameters
double kappa, tau, mu;
double par0, dpar;
double par[nparmax];

// observables per measurement
double energy[nmeasmax],mag[nmeasmax],nn[nmeasmax];

// linear contribution to the second derivatives of lnZ
double slin[nmeasmax],clin[nmeasmax], nlin[nmeasmax];
 
// average value of observables
double eaver[nparmax],eerr[nparmax]; // energy
double maver[nparmax],merr[nparmax]; // magnetization
double caver[nparmax],ccerr[nparmax]; // heat capacity
double saver[nparmax],serr[nparmax]; // susceptibility of mag.
double naver[nparmax],nerr[nparmax]; // particle number
double snaver[nparmax],snerr[nparmax]; // suscp. of particle number


// Subroutines 
void processdata( int &argc, char *argv[] );
void finalize(int it);
void calcdimerweights( );
void calcmonoweights( );

char texthelp[]="Usage: exec -f [FILE]\n"
		"Analysis program of the Z(3) spin model in the dual rep.\n"
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
  for (int ipar = 0; ipar<=npar ; ipar++ ){
  
  file << par[ipar]   << " "
       << eaver[ipar] << " " << eerr[ipar]  << " "
       << caver[ipar] << " " << ccerr[ipar] << " " 
       << maver[ipar] << " " << merr[ipar]  << " "
       << saver[ipar] << " " << serr[ipar]  << " " 
       << naver[ipar] << " " << nerr[ipar]  << " "
       << snaver[ipar] << " " << snerr[ipar] << " "
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
  int im;
  double esum,efsum,efblock,cfsum;
  double msum,mfsum,mfblock,sfsum;
  double nsum,nfsum,nfblock,snfsum;

  esum = 0.0;
  msum = 0.0;
  nsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    esum = esum + energy[im];
    msum = msum + mag[im];
    nsum = nsum + nn[im];
  }
  eaver[it] = esum/(double)nmeas;
  maver[it] = msum/(double)nmeas;
  naver[it] = nsum/(double)nmeas;

  efsum = 0.0;
  mfsum = 0.0;
  nfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    efsum += ( eaver[it] - energy[im] )*( eaver[it] - energy[im] );
    mfsum += ( maver[it] - mag[im] )*( maver[it] - mag[im] );
    nfsum += ( naver[it] - nn[im]  )*( naver[it] - nn[im]  );
  }
  eerr[it] = sqrt(efsum)/(double)nmeas;
  merr[it] = sqrt(mfsum)/(double)nmeas;
  nerr[it] = sqrt(nfsum)/(double)nmeas;
    
  for(im=0 ; im<nmeas ; im++){
    efsum += clin[im];
    mfsum += slin[im];
    nfsum += nlin[im];
  }
  caver[it] = efsum/(double)nmeas;
  saver[it] = mfsum/(double)nmeas;
  snaver[it] = nfsum/(double)nmeas;
    
  cfsum = 0.0;
  sfsum = 0.0;
  snfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock = ( efsum - ( eaver[it] - energy[im] )*( eaver[it] - energy[im] ) - clin[im] )
							/double(nmeas-1);
    mfblock = ( mfsum - ( maver[it] - mag[im]    )*( maver[it] - mag[im]    ) - slin[im] )
							/double(nmeas-1);
    nfblock = ( nfsum - ( naver[it] - nn[im] )*( naver[it] - nn[im] ) - nlin[im] )
							/double(nmeas-1);

    cfsum  += ( caver[it] - efblock )*( caver[it] - efblock );
    sfsum  += ( saver[it] - mfblock )*( saver[it] - mfblock );
    snfsum += ( snaver[it] - nfblock )*( snaver[it] - nfblock );
  }  
  ccerr[it] = sqrt(cfsum);
  serr[it] = sqrt(sfsum);
  snerr[it] = sqrt(snfsum);
   
  eaver[it] = eaver[it]/(double)nsite;
  eerr[it]  = eerr[it]/(double)nsite;
  caver[it] = caver[it]/(double)nsite;
  ccerr[it]  = ccerr[it]/(double)nsite;  

  maver[it] = maver[it]/(double)nsite;
  merr[it]  = merr[it]/(double)nsite;
  saver[it] = saver[it]/(double)nsite;
  serr[it]  = serr[it]/(double)nsite;  

  naver[it] = naver[it]/(double)nsite;
  nerr[it]  = nerr[it]/(double)nsite;
  snaver[it] = snaver[it]/(double)nsite;
  snerr[it]  = snerr[it]/(double)nsite;
}

//_________________________________________________________________________
inline void getfilename(  int &argc, char *argv[] )
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
   
//_________________________________________________________________________
void processdata( int &argc, char *argv[] )
{
	// local MC parameters
  int nequi,nskip,iseed;

	// link occupation number (dimers)
  int ndimp,ndimm,ndim0;

	// site occupation number (monomers)
	int nmonp,nmonm,nmon0;

	// weights that depend on the monomer number
  double rrm=0.,rr0=0.,rrp=0.,rrrrm=0.,rrrr0=0.,rrrrp=0.;
  double rm=0.,r0=0.,rp=0.;
  double ffm=0.,ff0=0.,ffp=0.,qqm=0.,qq0=0.,qqp=0.;
  double pm=0.,p0=0.,pp=0.,pbm=0.,pb0=0.,pbp=0.;
  double qm=0.,q0=0.,qp=0.,qbm=0.,qb0=0.,qbp=0.;

  int *nblock;
  
	getfilename( argc, argv );
  
  file.open(infile, ios::in | ios::binary);
  file.read((char*)&leng,sizeof(int));
#ifdef KAPPA  
  file.read((char*)&tau,sizeof(double));
  file.read((char*)&mu,sizeof(double));
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
#endif
#ifdef MU
  file.read((char*)&tau,sizeof(double));
  file.read((char*)&kappa,sizeof(double));
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
#endif
#ifdef TAU
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
  file.read((char*)&kappa,sizeof(double));
  file.read((char*)&mu,sizeof(double));
#endif
  file.read((char*)&nequi,sizeof(int));
  file.read((char*)&nmeas,sizeof(int));
  file.read((char*)&nskip,sizeof(int));
  file.read((char*)&iseed,sizeof(int));
   
  printf(" leng    = %d\n", leng);
#ifdef KAPPA
  printf(" kappa0  = %f\n", par0);
  printf(" dkappa  = %f\n", dpar);
  printf(" nkappa  = %d\n", npar);
  printf(" mu      = %f\n", mu);
  printf(" tau     = %f\n", tau);
#endif
#ifdef MU
  printf(" mu0     = %f\n", par0);
  printf(" dmu     = %f\n", dpar);
  printf(" nmu     = %d\n", npar);
  printf(" kappa   = %f\n", kappa);
  printf(" tau     = %f\n", tau);
#endif
#ifdef TAU
  printf(" tau0    = %f\n", par0);
  printf(" dtau    = %f\n", dpar);
  printf(" ntau    = %d\n", npar);
  printf(" kappa   = %f\n", kappa);
  printf(" mu      = %f\n", mu);
#endif
  printf(" nequi   = %d\n", nequi);
  printf(" nmeas   = %d\n", nmeas);
  printf(" nskip   = %d\n", nskip);
  printf(" iseed   = %d\n", iseed);

  nsite = leng*leng*leng;

  nblock = new int[6*nmeas];
  		  		   
  for(int ipar=0 ; ipar<=npar ; ipar++ )
	{  
    par[ipar] = par0 + ipar*dpar;
#ifdef KAPPA
		kappa = par[ipar];
#endif
#ifdef MU
		mu = par[ipar];
#endif
#ifdef TAU
		tau = par[ipar];
#endif

    double eta = kappa*exp(mu);
    double etab = kappa*exp(-mu);
		calcdimerweights( );
    calcmonoweights( );

    if (kappa > 0.0000000001){
    rrp = kappa*( exp(-mu)*monoweight[0] + exp(mu)*monoweight[1] )/monoweight[2];
    rr0 = kappa*( exp(-mu)*monoweight[2] + exp(mu)*monoweight[0] )/monoweight[1];
    rrm = kappa*( exp(-mu)*monoweight[1] + exp(mu)*monoweight[2] )/monoweight[0];

    rp = kappa*( -exp(-mu)*monoweight[0] + exp(mu)*monoweight[1] )/monoweight[2];
    r0 = kappa*( -exp(-mu)*monoweight[2] + exp(mu)*monoweight[0] )/monoweight[1];
    rm = kappa*( -exp(-mu)*monoweight[1] + exp(mu)*monoweight[2] )/monoweight[0];
		
    rrrrp = kappa*kappa*( 2*( monoweight[2]*monoweight[2] - monoweight[0]*monoweight[1] )
          + exp(-2*mu)*( monoweight[1]*monoweight[2] - monoweight[0]*monoweight[0] ) 
          + exp( 2*mu)*( monoweight[0]*monoweight[2] - monoweight[1]*monoweight[1] ) ) / 
	  ( monoweight[2]*monoweight[2] ) ;
          	
    rrrr0 = kappa*kappa*( 2*( monoweight[1]*monoweight[1] - monoweight[0]*monoweight[2] )
          + exp(-2*mu)*( monoweight[0]*monoweight[1] - monoweight[2]*monoweight[2] )
          + exp( 2*mu)*( monoweight[2]*monoweight[1] - monoweight[0]*monoweight[0] ) ) / 
	  ( monoweight[1]*monoweight[1] ) ; 
          	
    rrrrm = kappa*kappa*( 2*( monoweight[0]*monoweight[0] - monoweight[2]*monoweight[1] )
          + exp(-2*mu)*( monoweight[2]*monoweight[0] - monoweight[1]*monoweight[1] )
          + exp( 2*mu)*( monoweight[1]*monoweight[0] - monoweight[2]*monoweight[2] ) ) /
          ( monoweight[0]*monoweight[0] ) ;

    ffp = monoweight[1]/monoweight[2];
    ff0 = monoweight[0]/monoweight[1];
    ffm = monoweight[2]/monoweight[0];
  
    qqp = (monoweight[0]*monoweight[2] - monoweight[1]*monoweight[1])/(monoweight[2]*monoweight[2]);
    qq0 = (monoweight[2]*monoweight[1] - monoweight[0]*monoweight[0])/(monoweight[1]*monoweight[1]);
    qqm = (monoweight[1]*monoweight[0] - monoweight[2]*monoweight[2])/(monoweight[0]*monoweight[0]) ;

    pp = monoweight[1]/monoweight[2];
    p0 = monoweight[0]/monoweight[1];
    pm = monoweight[2]/monoweight[0];

    pbp = monoweight[0]/monoweight[2];
    pb0 = monoweight[2]/monoweight[1];
    pbm = monoweight[1]/monoweight[0];

    rp = 1. - pp*pbp;
    r0 = 1. - p0*pb0;
    rm = 1. - pm*pbm;
    // for F(mu)
    //rp = kappa*( -exp(-mu)*monoweight[0] + exp(mu)*monoweight[1] )/monoweight[2];
    //r0 = kappa*( -exp(-mu)*monoweight[2] + exp(mu)*monoweight[0] )/monoweight[1];
    //rm = kappa*( -exp(-mu)*monoweight[1] + exp(mu)*monoweight[2] )/monoweight[0];

    qbp = monoweight[1]/monoweight[2] - pbp*pbp;
    qb0 = monoweight[0]/monoweight[1] - pb0*pb0;
    qbm = monoweight[2]/monoweight[0] - pbm*pbm;

    qp = monoweight[0]/monoweight[2] - pp*pp;
    q0 = monoweight[2]/monoweight[1] - p0*p0;
    qm = monoweight[1]/monoweight[0] - pm*pm;
    }

    file.read( (char*)nblock, 6*nmeas*sizeof(int) );

    for ( int imeas = 0; imeas<nmeas; imeas++ )
    {
       ndimm = nblock[imeas*6+0];
       ndim0 = nblock[imeas*6+1];
       ndimp = nblock[imeas*6+2];
       
       nmonm = nblock[imeas*6+3];
       nmon0 = nblock[imeas*6+4];
       nmonp = nblock[imeas*6+5];

       energy[imeas] = -3.0*nsite*cc - bb*(ndimm + ndimp) - (nmonm*rrm + nmon0*rr0 + nmonp*rrp);
       mag[imeas] = nmonm*ffm + nmon0*ff0 + nmonp*ffp;
       double magb = nmonm*pbm + nmon0*pb0 + nmonp*pbp;
       nn[imeas] = eta*mag[imeas] - etab*magb;

       clin[imeas] = 3.0*nsite*cccc + bbbb*(ndimm + ndimp) + nmonm*rrrrm + nmon0*rrrr0 + nmonp*rrrrp;
       slin[imeas] = nmonm*qqm + nmon0*qq0 + nmonp*qqp;
       double slinb = nmonm*qbm + nmon0*qb0 + nmonp*qbp;
       nlin[imeas] = eta*mag[imeas] + etab*magb +
                     eta*eta*slin[imeas] + etab*etab*slinb +
                     -2.*eta*etab*( nmonm*rm + nmon0*r0 + nmonp*rp );
    }
    
    cout << "reading tau=" << tau << " - kappa=" << kappa << " - mu=" << mu << endl;
    
    finalize(ipar);
    
  }
  
  file.close();
}

//________________________________________________________________________
void calcdimerweights()
{
  cc = 2.*tau*( exp(2.*tau)-exp(-tau) )/( exp(2.*tau) + 2.*exp(-tau) );
  cccc = 18.*tau*tau*exp(tau) / ( ( exp(2.*tau) + 2.*exp(-tau) )*( exp(2.*tau) + 2.*exp(-tau) ) );
  bb = 9.*tau/( exp(3.*tau) + 1. - 2.*exp(-3.*tau) );
  bbbb = -tau*tau*( 27.*exp(3.*tau) + 54.*exp(-3.*tau) ) 
				/( ( exp(3.*tau) + 1. - 2.*exp(-3.*tau) )*( exp(3.*tau) + 1. - 2.*exp(-3.*tau) ) );
}

//_________________________________________________________________________
void calcmonoweights( )
{
  int l;
  const double pi2 = 2.0*3.141592653589793;

  monoweight[0] = 0.;
  monoweight[1] = 0.;
  monoweight[2] = 0.;
    
  for (l = -1; l<=1; ++l){  
    monoweight[l+1] = ( exp(2.0*kappa*cosh(mu)) + 2.0*exp(-kappa*cosh(mu))*
                      cos(sqrt(3.0)*kappa*sinh(mu) - l*pi2/3.0) )/3.0;
  }
}
