/************************************************************************
* Analysis program for different values of mu 
* in the dual representation of the SU(3) spin model (new variables).  
* See ref.: arXiv:1204.6074
*
* To execute: ./bin/anal_rot_tau.x -f file_with_configs
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <getopt.h>
using namespace std;

fstream file;
char infile[150],outfile[150];

// MAX. dimension of arrays
const int nparmax  = 41;
const int nmeasmax = 1000000;

// Lattice size
int leng,nsite;

// Global MC parameters
int npar,nmeas;
double tau,mu,kappa,eta,etabar;
double par[nparmax];
double par0,dpar;

// observables per measurement
double energy[nmeasmax],mag[nmeasmax];

// linear contribution to the second derivatives of lnZ
double clin[nmeasmax],slin[nmeasmax];
 
// average value of observables
double eaver[nparmax],eerr[nparmax]; // energy
double maver[nparmax],merr[nparmax]; // magnetization
double caver[nparmax],ccerr[nparmax]; // heat capacity
double saver[nparmax],serr[nparmax]; // susc. of mag.


// Subroutines
void processdata( int &argc, char *argv[] );
void finalize(int it);

char texthelp[]="Usage: exec -f [FILE]\n"
    "Analysis program of the SU(3) spin model in the dual rep. (new var.)\n"
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
    file << par[ipar]    << " " 
         << eaver[ipar]  << " " << eerr[ipar]  << " "
         << caver[ipar]  << " " << ccerr[ipar] << " " 
         << maver[ipar]  << " " << merr[ipar]  << " "
         << saver[ipar]  << " " << serr[ipar] 
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

  esum = 0.0;
  msum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    esum = esum + energy[im];
    msum = msum + mag[im];
  }
  eaver[it] = esum/(double)nmeas;
  maver[it] = msum/(double)nmeas;

  efsum = 0.0;
  mfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    efsum += ( eaver[it] - energy[im] )*( eaver[it] - energy[im] );
    mfsum += ( maver[it] - mag[im] )*( maver[it] - mag[im] );
  }
  eerr[it] = sqrt(efsum)/(double)nmeas;
  merr[it] = sqrt(mfsum)/(double)nmeas;

  for(im=0 ; im<nmeas ; im++){
    efsum += clin[im];
    mfsum += slin[im];
  }
  caver[it]  = efsum/(double)nmeas;    
  saver[it]  = mfsum/(double)nmeas;

  cfsum  = 0.0;
  sfsum  = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock = ( efsum - ( eaver[it] - energy[im] )*( eaver[it] - energy[im] ) - clin[im] )/double(nmeas-1);
    mfblock = ( mfsum - ( maver[it] - mag[im] )*( maver[it] - mag[im] ) - slin[im] )/double(nmeas-1);

    cfsum  += ( caver[it] - efblock )*( caver[it] - efblock );      
    sfsum  += ( saver[it] - mfblock )*( saver[it] - mfblock );
  }  
  ccerr[it] = sqrt(cfsum);
  serr[it]  = sqrt(sfsum);
   
  eaver[it] = eaver[it]/(double)nsite;
  eerr[it]  = eerr[it]/(double)nsite;
  caver[it] = caver[it]/(double)nsite;
  ccerr[it] = ccerr[it]/(double)nsite;  

  maver[it] = maver[it]/(double)nsite;
  merr[it]  = merr[it]/(double)nsite;
  saver[it] = saver[it]/(double)nsite;
  serr[it]  = serr[it]/(double)nsite;
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

  // Monomer and dimer occupation numbers
  int nmon[3], ndim;

  // Array with measurements from gen.cpp
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
   
  printf(" leng  = %d\n", leng);
#ifdef MU
  printf(" mu0   = %f\n", par0);
  printf(" dmu   = %f\n", dpar);
  printf(" nmu   = %d\n", npar);
  printf(" kappa = %f\n", kappa);
  printf(" tau   = %f\n", tau);
#else
  printf(" tau0   = %f\n", par0);
  printf(" dtau   = %f\n", dpar);
  printf(" ntau   = %d\n", npar);
  printf(" kappa  = %f\n", kappa);
  printf(" mu     = %f\n", mu);
#endif
  printf(" nequi  = %d\n", nequi);
  printf(" nmeas  = %d\n", nmeas);
  printf(" nskip  = %d\n", nskip);
  printf(" iseed  = %d\n", iseed);

  nsite  = leng*leng*leng;
  nblock = new int[4*nmeas];
   
  for(int ipar=0 ; ipar<=npar ; ipar++ )
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
    eta = 1.0/(kappa*exp(mu));
    etabar = 1.0/(kappa*exp(-mu));
#endif

    file.read( (char*)nblock, 4*nmeas*sizeof(int) );

    for ( int imeas = 0; imeas<nmeas; imeas++ ){
       ndim = nblock[imeas*4+0];       
       nmon[0] = nblock[imeas*4+1];
       nmon[1] = nblock[imeas*4+2];
       nmon[2] = nblock[imeas*4+3];
   
       energy[imeas] = ndim + nmon[2] + nmon[0];
       clin[imeas]   = -energy[imeas];

       mag[imeas]  = 0.5*(nmon[0] + nmon[1] + nmon[2])*eta;
       slin[imeas] = -mag[imeas]*eta;
    }
    cout << "Reading tau= " << tau << " - kappa= " << kappa << " - mu= " << mu << endl; 

    finalize(ipar);  
  }
  
  file.close();
  
  //delete[] nblock;
   
}  
