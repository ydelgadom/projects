/************************************************************************
* Analysis program for different values of beta 
*	in the dual representation of the U(1) gauge-Higgs model (2 flavors).  
*	See references in ../papers.
*
*	To execute: ./bin/anal.x -f file_with_configs (without id number and without extension)
* eg. if config files are:
*			file_0.out
*			file_1.out
* --> execute:  ./bin/anal.x -f PATH_TO_FILES/file 
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

#define ANALYSIS
#define LENGTH 50
#include "weights.h"

char infile[128],outfile[128];

const int nbetamax = 40;
const int nmeasmax = 500000;

int leng,leng_t,nsite;
int nmeas0,nmeas[nbetamax];
int nvar;  // number of configuration files
int ncut;	 // to discard first ncut measurements

double kappa0,dkappa,lambda,mu0,dmu,beta0,dbeta;

// Observables per measurment
double u[nmeasmax],phi2[2][nmeasmax],nlink[2][nmeasmax];

// Linear contributions to the second derivative
double clin[nmeasmax],sphi2lin[2][nmeasmax];
  
// Average value of observables
double uaver[nbetamax],uerr[nbetamax]; // derivative w.r.t. beta
double caver[nbetamax],ccerr[nbetamax];	// 2nd derivative w.r.t. beta
double phi2aver[2][nbetamax],phi2err[2][nbetamax];	// derivative w.r.t kappa
double sphi2aver[2][nbetamax],sphi2err[2][nbetamax];	// 2nd. derivative w.r.t kappa
double linkaver[2][nbetamax],linkerr[2][nbetamax];	// link occupation number
double slinkaver[2][nbetamax],slinkerr[2][nbetamax];	// susceptibility of link occup. number

fstream file;

void processdata( int &argc, char *argv[]  );
void finalize( int it );

char texthelp[]="Usage: exec -f [FILE]\n"
		"Analysis program of the U(1) gauge-Higgs model in the dual rep.\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -f, --F  File with the configurations\n"
		"Report bugs to ydelgado83@gmail.com\n";

//-----------------------------------------------------------------------------   
int main( int argc, char *argv[] )
{
  printf( "\nProgram analysis.cpp\n" );

  processdata( argc, argv );

	// print observables in output file   
  file.open( outfile, ios::out | ios::trunc );
  int i = 0;
  double betan = beta0 + dbeta*nvar;
  double mun = mu0 + dmu*nvar;
  double kappan = kappa0 + dkappa*nvar;
  for ( double beta=beta0 ; beta<=betan ; beta+=dbeta )
  {
    for ( double kappa=kappa0 ; kappa<=kappan ; kappa+=dkappa )
    {
      for ( double mu=mu0 ; mu<=mun ; mu+=dmu )
      {
				if (i==nvar) break;

        file << beta << " " << kappa << " " << mu << " "
             << uaver[i]      << " " << uerr[i]     << " "
             << caver[i]      << " " << ccerr[i]    << " "

             << phi2aver[0][i]   << " " << phi2err[0][i]  << " "
             << sphi2aver[0][i]  << " " << sphi2err[0][i] << " "
             << phi2aver[1][i]   << " " << phi2err[1][i]  << " "
             << sphi2aver[1][i]  << " " << sphi2err[1][i] << " "

             << linkaver[0][i]   << " " << linkerr[0][i]  << " "
             << slinkaver[0][i]  << " " << slinkerr[0][i] << " "
             << linkaver[1][i]   << " " << linkerr[1][i]  << " "
             << slinkaver[1][i]  << " " << slinkerr[1][i] << " "
             << endl;
      
        ++i;
        if ( dmu<=0 ) break;
      }
      if ( dkappa<=0 ) break;
    }
    if ( dbeta <= 0 ) break;
  }
  file.close( );
  printf( "\nDone.\n" );
}
  

//!------------------------------------------------------------------------
//!------------------------------------------------------------------------
void finalize( int it )
{
	/*
		This subrouting computes the average values of the
		observables and the errors (Jackknife is used
		for the susceptibilities)
	*/ 
  int im;
  double efblock, lfblock[2], mfblock[2];

  double esum = 0.0;
  double msum[2] = {0.0,0.0};
  double lsum[2] = {0.0,0.0};
  for(im=ncut ; im<nmeas[it] ; im++)
  {
    esum = esum + u[im];

    msum[0] += phi2[0][im];
    lsum[0] += nlink[0][im];

    msum[1] += phi2[1][im];
    lsum[1] += nlink[1][im];
  }
  uaver[it] = esum/double(nmeas[it]-ncut);
  phi2aver[0][it] = msum[0]/double(nmeas[it]-ncut);
  linkaver[0][it] = lsum[0]/double(nmeas[it]-ncut);
  phi2aver[1][it] = msum[1]/double(nmeas[it]-ncut);
  linkaver[1][it] = lsum[1]/double(nmeas[it]-ncut);

  double efsum = 0.0;
  double mfsum[2] = {0.0,0.0};
  double lfsum[2] = {0.0,0.0};
  for( im=ncut ; im<nmeas[it] ; im++ )
  {
    efsum += pow( uaver[it] - u[im] , 2 );
    mfsum[0] += pow( phi2aver[0][it] - phi2[0][im] , 2 );
    lfsum[0] += pow( linkaver[0][it] - nlink[0][im] , 2 );

    mfsum[1] += pow( phi2aver[1][it] - phi2[1][im] , 2 );
    lfsum[1] += pow( linkaver[1][it] - nlink[1][im] , 2 );
  }
  uerr[it] = sqrt(efsum)/double(nmeas[it]-ncut);
  phi2err[0][it] = sqrt(mfsum[0])/double(nmeas[it]-ncut);
  linkerr[0][it] = sqrt(lfsum[0])/double(nmeas[it]-ncut);
  phi2err[1][it] = sqrt(mfsum[1])/double(nmeas[it]-ncut);
  linkerr[1][it] = sqrt(lfsum[1])/double(nmeas[it]-ncut);

  for( im=ncut ; im<nmeas[it] ; im++ )
  {
    efsum += clin[im];
    mfsum[0] += sphi2lin[0][im];
    mfsum[1] += sphi2lin[1][im];
  }
  caver[it] = efsum/double(nmeas[it]-ncut);
  sphi2aver[0][it] = mfsum[0]/double(nmeas[it]-ncut);
  slinkaver[0][it] = lfsum[0]/double(nmeas[it]-ncut);
  sphi2aver[1][it] = mfsum[1]/double(nmeas[it]-ncut);
  slinkaver[1][it] = lfsum[1]/double(nmeas[it]-ncut);

  double cfsum = 0.0;
  double sfsum[2] = {0.0,0.0};
  double slsum[2] = {0.0,0.0};
  for( im=ncut ; im<nmeas[it] ; im++ )
  {  
    efblock  = ( efsum - pow( uaver[it] - u[im], 2 ) - clin[im] )
								/double(nmeas[it]-ncut-1);
    mfblock[0]  = ( mfsum[0] - pow( phi2aver[0][it] - phi2[0][im], 2 ) - sphi2lin[0][im] )
								/double(nmeas[it]-ncut-1);
    lfblock[0]  = ( lfsum[0] - pow( linkaver[0][it] - nlink[0][im], 2 ) )
								/double(nmeas[it]-ncut-1);
    mfblock[1]  = ( mfsum[1] - pow( phi2aver[1][it] - phi2[1][im], 2 ) - sphi2lin[1][im] )
								/double(nmeas[it]-ncut-1);
    lfblock[1]  = ( lfsum[1] - pow( linkaver[1][it] - nlink[1][im], 2 ) )
								/double(nmeas[it]-ncut-1);

    cfsum  += pow( caver[it] - efblock , 2 );
    sfsum[0]  += pow( sphi2aver[0][it] - mfblock[0] , 2 );
    slsum[0]  += pow( slinkaver[0][it] - lfblock[0] , 2 );
    sfsum[1]  += pow( sphi2aver[1][it] - mfblock[1] , 2 );
    slsum[1]  += pow( slinkaver[1][it] - lfblock[1] , 2 );
  }  
  ccerr[it] = sqrt(cfsum);
  sphi2err[0][it] = sqrt(sfsum[0]);
  slinkerr[0][it] = sqrt(slsum[0]);
  sphi2err[1][it] = sqrt(sfsum[1]);
  slinkerr[1][it] = sqrt(slsum[1]);
   
  phi2aver[0][it]  /= (double)nsite;
  phi2err[0][it]   /= (double)nsite;
  sphi2aver[0][it] /= (double)nsite;
  sphi2err[0][it]  /= (double)nsite;

  phi2aver[1][it]  /= (double)nsite;
  phi2err[1][it]   /= (double)nsite;
  sphi2aver[1][it] /= (double)nsite;
  sphi2err[1][it]  /= (double)nsite;

  uaver[it]  /= (double)(6*nsite);
  uerr[it]   /= (double)(6*nsite);
  caver[it]  /= (double)(6*nsite);
  ccerr[it]  /= (double)(6*nsite);

  linkaver[0][it]  /= (double)(nsite);
  linkerr[0][it]   /= (double)(nsite);
  slinkaver[0][it] /= (double)(nsite);
  slinkerr[0][it]  /= (double)(nsite);

  linkaver[1][it]  /= (double)(nsite);
  linkerr[1][it]   /= (double)(nsite);
  slinkaver[1][it] /= (double)(nsite);
  slinkerr[1][it]  /= (double)(nsite);

}
   
//!------------------------------------------------------------------------
void processdata( int &argc, char *argv[] )
{  
  int nequi,nskip,sleng;
  char filename[100];
  
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
				sprintf(filename,"%s_0.out",optarg);
			  sprintf(outfile,"%s.obs",infile);
				cout << "Input File: " << infile << endl;
  			cout << "Output file: " << outfile << endl;
				break;

			default:
				cout << endl << texthelp << endl;
		}
	}
  file.open( filename, ios::in | ios::binary );

  file.read((char*)&leng,sizeof(int));
  file.read((char*)&sleng,sizeof(int));
  file.read((char*)&leng_t,sizeof(int));
  file.read((char*)&nvar,sizeof(int));

  file.read((char*)&lambda,sizeof(double));

  file.read((char*)&beta0,sizeof(double));
  file.read((char*)&dbeta,sizeof(double));

  file.read((char*)&kappa0,sizeof(double));
  file.read((char*)&dkappa,sizeof(double));

  file.read((char*)&mu0,sizeof(double));
  file.read((char*)&dmu,sizeof(double));

  file.read((char*)&nequi,sizeof(int));
  file.read((char*)&nmeas0,sizeof(int));
  file.read((char*)&nskip,sizeof(int));
     
  printf(" leng   = %d\n", leng);
  printf(" leng   = %d\n", sleng);
  printf(" leng_t = %d\n", leng_t);
  printf(" nvar   = %d\n", nvar );
  printf(" lambda = %f\n", lambda);

  printf(" beta0  = %f\n", beta0);
  printf(" dbeta  = %f\n", dbeta);

  printf(" kappa0 = %f\n", kappa0);
  printf(" dkappa = %f\n", dkappa);

  printf(" mu0    = %f\n", mu0);
  printf(" dmu    = %f\n", dmu);

  printf(" nequi  = %d\n", nequi);
  printf(" nmeas  = %d\n", nmeas0);
  printf(" nskip  = %d\n", nskip);
  
  nsite  = leng*leng*leng*leng_t;

  int fcc = leng/sleng;
  int snsite = sleng*sleng*sleng*leng_t;
  fcc = fcc*fcc*fcc;

  cout << "\nInsert cut for measurements:\n" << endl;
  cin >> ncut;

  int nblock[500], nbytes;
  for ( int ivar=0; ivar<nvar; ivar++ )
  {
    double beta = beta0 + dbeta*ivar;
    double mu = mu0 + dmu*ivar;
    double kappa = kappa0 + dkappa*ivar;

    calculate_Ibeta( beta );
    calculate_Pn( lambda, kappa );

    sprintf( filename, "%s_%d.out", infile, ivar);
    cout << "reading ... " << filename <<  endl;
    if ( ivar!=0 ) file.open( filename, ios::in | ios::binary );
        
    nmeas[ivar] = nmeas0;
    for ( int imeas = 0; imeas<nmeas0; imeas++ )
    {
      file.read( (char*)&nbytes, sizeof(int) );
      file.read( (char*)nblock, nbytes*sizeof(int) );

      int inp=0, inf[2]={0,0};
      int is_plaq=1, is_flux0=0, next=1;
      int total_np=0, total_nf[2]={0,0};
      int ib;

      u[imeas] = 0.0;
      clin[imeas] = 0.0;
      phi2[0][imeas] = 0.0; phi2[1][imeas] = 0.0;
      sphi2lin[0][imeas] = 0.0; sphi2lin[1][imeas] = 0.0;

      nlink[0][imeas] = nblock[0];
      nlink[1][imeas] = nblock[1];

      for ( ib=nbytes; ib>=0; ib-- ){
        if(is_plaq){
          int np = nblock[2+inp];
          if ( np != -1 ){
            if (inp==0){
              np -= (nsite - snsite)*6*fcc;
            }
            u[imeas] += np * dIn[inp];
            clin[imeas] += np * ( d2In[inp] - pow(dIn[inp], 2) );
            total_np += np;
            inp++;
          }
          else{
            is_plaq=0;
            is_flux0=1;
           }
        }
        else if (is_flux0){
         int nf = nblock[inp+3+inf[0]];
         if ( nf != -2 ){
            if ( inf[0]==0 ){
              nf -= (nsite-snsite)*fcc;
            }
            phi2[0][imeas] += nf * ( Pn[inf[0]+2]/Pn[inf[0]] );
            sphi2lin[0][imeas] += nf * ( Pn[inf[0]+4]/Pn[inf[0]] - pow( Pn[inf[0]+2]/Pn[inf[0]], 2 ) );
            total_nf[0] += nf;
            inf[0]++;
          }
          else{
            is_flux0=0;
          }
        }
        else{
          int nf = nblock[inp+4+inf[0]+inf[1]];
          if ( nf != -3 ){
            if ( inf[1]==0 ){
              nf -= (nsite-snsite)*fcc;
            }
            phi2[1][imeas] += nf  * ( Pn[inf[1]+2]/Pn[inf[1]] );
            sphi2lin[1][imeas] += nf  * ( Pn[inf[1]+4]/Pn[inf[1]] - pow( Pn[inf[1]+2]/Pn[inf[1]], 2 ) );
            total_nf[1] += nf;
            inf[1]++;
          }
          else{
            next=0;
          }
        }
        if (next==0) break;
      }
      /*if (nbytes!=(inp+inf[0]+inf[1]+5)){
        cout << "ERROR reading data " << imeas << endl;
        nmeas[ivar] = imeas; break;
      }*/
      if ( total_np != 6*nsite ){
        cout << "ERROR: np != 6*nsite " << total_np <<"/"<< 6*nsite << " " << imeas << endl;
        nmeas[ivar] = imeas; break;
      }
      if ( total_nf[0] != nsite ){
        cout << "ERROR: nf0 != nsite " << total_nf[0] << "/" << nsite << " " << imeas << endl;
        nmeas[ivar] = imeas; break;
      }
      if ( total_nf[1] != nsite ){
        cout << "ERROR: nf1 != nsite " << total_nf[1] << "/" << nsite << " " << imeas << endl;
        nmeas[ivar] = imeas; break;
      }
      if ( file.eof( ) ) { 
        nmeas[ivar] = imeas; 
        cout << "end of file " << imeas << endl; 
        break; 
      }
     }
     printf("%d: (b, k, m) = (%f, %f, %f) - nmeas = %d\n", ivar, beta, kappa, mu, nmeas[ivar]);
     finalize( ivar );
     file.close( );
  }
}
