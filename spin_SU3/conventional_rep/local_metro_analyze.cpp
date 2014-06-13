/************************************************************************
* Analysis program for the conventional representation of the 
*	SU(3) spin model.  See ref.: arXiv:1204.6074
*
*	To execute: ./bin/anal.x -f file_with_configs
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

const int ntaumax  = 41;
const int nmeasmax = 1000000;

int    leng,nsite,ntau,nmeas;
double eps,kappa;
double tau[ntaumax];

// observables per measurement
double energy[nmeasmax],mag[nmeasmax];
 
// average value of observables
double eaver[ntaumax],eerr[ntaumax]; // energy
double maver[ntaumax],merr[ntaumax]; // magnetization
double caver[ntaumax],ccerr[ntaumax]; // heat capacity
double saver[ntaumax],serr[ntaumax]; // susceptibility of mag.

fstream file;

void processdata( int &argc, char *argv[] );
void finalize(int it);

char texthelp[]="Usage: exec -f [FILE]\n"
		"Analysis program of the SU(3) spin model in the conventional rep.\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -f, --F  File with the configurations\n"
		"Report bugs to ydelgado83@gmail.com\n";

//-----------------------------------------------------------------------------   
int main( int argc, char *argv[] )
{
  printf("\nProgram metro_su3_analyze.c\n");

  processdata( argc, argv );

	// print observables in output file
  file.open(outfile, ios::out | ios::trunc );
  for (int itau = 0; itau<=ntau ; itau++ )
  {
    file << tau[itau]    << " "
         << eaver[itau]  << " " << eerr[itau]  << " "
         << caver[itau]  << " " << ccerr[itau] << " " 
         << maver[itau]  << " " << merr[itau]  << " "
         << saver[itau]  << " " << serr[itau] 
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

  caver[it]  = efsum/(double)nmeas;    
  saver[it]  = mfsum/(double)nmeas;

  cfsum  = 0.0;
  sfsum  = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock = ( efsum - ( eaver[it] - energy[im] )*( eaver[it] - energy[im] ) )/double(nmeas-1);
    mfblock = ( mfsum - ( maver[it] - mag[im] )*( maver[it] - mag[im] ) )/double(nmeas-1);

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
void processdata( int &argc, char *argv[] )
{  
  int    nequi,nskip,iseed,itau,imeas;
  double tau0,dtau;
  double *nblock;

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
  
  file.open(infile, ios::in | ios::binary);
  
  file.read((char*)&leng,sizeof(int));
  file.read((char*)&tau0,sizeof(double));
  file.read((char*)&dtau,sizeof(double));
  file.read((char*)&ntau,sizeof(int));
  file.read((char*)&kappa,sizeof(double));
  file.read((char*)&eps,sizeof(double));
  file.read((char*)&nequi,sizeof(int));
  file.read((char*)&nmeas,sizeof(int));
  file.read((char*)&nskip,sizeof(int));
  file.read((char*)&iseed,sizeof(int));
   
  printf(" leng    = %d\n", leng);
  printf(" tau0    = %f\n", tau0);
  printf(" dtau    = %f\n", dtau);
  printf(" ntau    = %d\n", ntau);
  printf(" kappa   = %f\n", kappa);
  printf(" eps     = %f\n", eps);
  printf(" nequi   = %d\n", nequi);
  printf(" nmeas   = %d\n", nmeas);
  printf(" nskip   = %d\n", nskip);
  printf(" iseed   = %d\n", iseed);
  
  nsite  = leng*leng*leng;

  nblock = new double[2*nmeas];
  		  		   
  for(itau=0 ; itau<=ntau ; itau++ )
  {    
    tau[itau] = tau0 + itau*dtau;

    file.read( (char*)nblock, 2*nmeas*sizeof(double) );

    for ( imeas = 0; imeas<nmeas; imeas++ )
    {
       energy[imeas] = nblock[imeas*2+0];       
       mag[imeas]    = nblock[imeas*2+1];
    }
    
    printf("reading tau = %f\n", tau[itau]);
    
    finalize(itau);
    
  }
  
  file.close();
   
}  
