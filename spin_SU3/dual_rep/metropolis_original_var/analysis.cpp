/************************************************************************
* Analysis program for different values of mu 
* in the dual representation of the SU(3) spin model (before changing variables).  
* See ref.: arXiv:1204.6074
*
* To execute: ./bin/anal_mu.x -f file_with_configs
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
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

// Array dimensions
const int nparmax   = 41;
const int nmeasmax = 1000000;

// Lattice size
int leng,nsite;

// Global MC parameters
int nmeas;

// Global parameters
double kappa,mu,tau,eta,etabar;
double par0,dpar;
double par[nparmax];
int npar;

// observables per measurement
double pp[nmeasmax];
double energy[nmeasmax],mag[nmeasmax];

// linear contribution to the second derivatives of lnZ
double clin[nmeasmax],slin[nmeasmax],plin[nmeasmax];
  
// average observables
double eaver[nparmax],eerr[nparmax];  // energy
double maver[nparmax],merr[nparmax];  // magnetization
double caver[nparmax],ccerr[nparmax]; // heat capacity
double saver[nparmax],serr[nparmax];  // susceptibility of mag.
double paver[nparmax],perr[nparmax];  
double spaver[nparmax],sperr[nparmax];


// Subroutines
void processdata( int &argc, char *argv[] );
void finalize(int it);

char texthelp[]="Usage: exec -f [FILE]\n"
    "Analysis program of the SU(3) spin model in the dual rep. (original var.)\n"
    "\n"
    "Mandatory arguments to long options are mandatory for short options too.\n"
    "  -f, --F  File with the configurations\n"
    "Report bugs to ydelgado83@gmail.com\n";

//-----------------------------------------------------------------------------   
int main( int argc, char *argv[] )
{
  printf("\nProgram analysis.c\n");

  processdata( argc, argv );

  // print data to outputfile
  file.open(outfile, ios::out | ios::trunc );
  for (int ipar = 0; ipar<=npar ; ipar++ ){
    file << par[ipar]    << " "
         << eaver[ipar]  << " " << eerr[ipar]  << " "
         << caver[ipar]  << " " << ccerr[ipar] << " "
         << maver[ipar]  << " " << merr[ipar]  << " "
         << saver[ipar]  << " " << serr[ipar]  << " "
         << paver[ipar]  << " " << perr[ipar]  << " "
         << spaver[ipar] << " " << sperr[ipar]
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
  double psum,pfsum,pfblock,spfsum;

  esum = 0.0;
  msum = 0.0;
  psum = 0.0;

  for(im=0 ; im<nmeas ; im++){
    esum = esum + energy[im];
    msum = msum + mag[im];
    psum = psum + pp[im];
  }
  eaver[it] = esum/(double)nmeas;
  maver[it] = msum/(double)nmeas;
  paver[it] = psum/(double)nmeas;

  efsum = 0.0;
  mfsum = 0.0;
  pfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    efsum += ( eaver[it] - energy[im] )*( eaver[it] - energy[im] );
    mfsum += ( maver[it] - mag[im] )*( maver[it] - mag[im] );
    pfsum += ( paver[it] - pp[im] )*( paver[it] - pp[im] );
  }
  eerr[it] = sqrt(efsum)/(double)nmeas;
  merr[it] = sqrt(mfsum)/(double)nmeas;
  perr[it] = sqrt(pfsum)/(double)nmeas;

  for(im=0 ; im<nmeas ; im++){
    efsum += clin[im];
    mfsum += slin[im];
    pfsum += plin[im];
  }
  caver[it]  = efsum/(double)nmeas;    
  saver[it]  = mfsum/(double)nmeas;
  spaver[it] = pfsum/(double)nmeas;

  spfsum = 0.0;    
  cfsum  = 0.0;
  sfsum  = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock = ( efsum - ( eaver[it] - energy[im] )*( eaver[it] - energy[im] )  - clin[im] )/double(nmeas-1);
    mfblock = ( mfsum - ( maver[it] - mag[im] )*( maver[it] - mag[im] ) - slin[im] )/double(nmeas-1);
    pfblock = ( pfsum - ( paver[it] - pp[im] )*( paver[it] - pp[im] ) - plin[im] )/double(nmeas-1);

    cfsum  += ( caver[it] - efblock )*( caver[it] - efblock );      
    sfsum  += ( saver[it] - mfblock )*( saver[it] - mfblock );
    spfsum += ( spaver[it] - pfblock )*( spaver[it] - pfblock );
  }  
  ccerr[it] = sqrt(cfsum);
  serr[it]  = sqrt(sfsum);
  sperr[it] = sqrt(spfsum);
   
  eaver[it] = eaver[it]/(double)nsite;
  eerr[it]  = eerr[it]/(double)nsite;
  caver[it] = caver[it]/(double)nsite;
  ccerr[it] = ccerr[it]/(double)nsite;  

  maver[it] = maver[it]/(double)nsite;
  merr[it]  = merr[it]/(double)nsite;
  saver[it] = saver[it]/(double)nsite;
  serr[it]  = serr[it]/(double)nsite;

  paver[it]  = paver[it]/(double)nsite;
  perr[it]   = perr[it]/(double)nsite;
  spaver[it] = spaver[it]/(double)nsite;
  sperr[it]  = sperr[it]/(double)nsite;
   
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
  int nmon[3],ndim;

  // array with measurements from gen.cpp 
  int *nblock;
  
  getfilename( argc, argv );
  
  file.open(infile, ios::in | ios::binary);
  file.read((char*)&leng,sizeof(int));
#ifdef MU
  file.read((char*)&par0,sizeof(double));
  file.read((char*)&dpar,sizeof(double));
  file.read((char*)&npar,sizeof(int));
  file.read((char*)&kappa,sizeof(double));
  file.read((char*)&tau,sizeof(double));
#else
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
   
  printf(" leng   = %d\n", leng);
#ifdef MU
  printf(" mu0    = %f\n", par0);
  printf(" dmu    = %f\n", dpar);
  printf(" nmu    = %d\n", npar);
  printf(" kappa  = %f\n", kappa);
  printf(" tau    = %f\n", tau);
#else
  printf(" tau0    = %f\n", par0);
  printf(" dtau    = %f\n", dpar);
  printf(" ntau    = %d\n", npar);
  printf(" kappa   = %f\n", kappa);
  printf(" mu      = %f\n", mu);
#endif
  printf(" nequi  = %d\n", nequi);
  printf(" nmeas  = %d\n", nmeas);
  printf(" nskip  = %d\n", nskip);
  printf(" iseed  = %d\n", iseed);
  
  nsite  = leng*leng*leng;
  nblock = new int[3*nmeas];
               
  for( int ipar=0 ; ipar<=npar ; ipar++ )
  {      
    par[ipar] = par0 + ipar*dpar;
#ifdef MU
    mu = par[ipar];
#else
    tau = par[ipar];
#endif
#ifdef KAPPA0
    eta = 0.;
    etabar = 0.;
#else
    eta = 1./(kappa*exp(mu));
    etabar = 1./(kappa*exp(-mu));
#endif
    file.read( (char*)nblock, 3*nmeas*sizeof(int) );

    for ( int imeas = 0; imeas<nmeas; imeas++ ){
       ndim = nblock[imeas*3+0];       
       nmon[0] = nblock[imeas*3+1];
       nmon[1] = nblock[imeas*3+2];
   
       energy[imeas] = ndim + nmon[0] + nmon[1];
       pp[imeas]     = 0.5*(nmon[1]*etabar + nmon[0]*eta);
       mag[imeas]    = nmon[0]*eta;

       clin[imeas] = -energy[imeas];
       plin[imeas] = -pp[imeas]*etabar;
       slin[imeas] = -mag[imeas]*eta;

    }
    cout << "Reading tau= " << tau << " - kappa= " << kappa << " - mu= " << mu << endl;
 
    finalize(ipar);  
  }
  
  file.close();
  
  //delete[] nblock;
}  
