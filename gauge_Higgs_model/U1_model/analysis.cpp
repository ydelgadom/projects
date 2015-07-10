/************************************************************************
* Analysis program for different values of beta 
* in the dual representation of the U(1) gauge-Higgs model (only 1 flavor).  
* See references in ../papers.
*
* To execute: ./bin/anal.x -f file_with_configs
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

#define ANALYSIS
#define LENGTH 50
#include "weights.h"

char infile[128],outfile[128];

const int nbetamax = 50;
const int nmeasmax = 1000000;

int    leng,leng_t,nsite,nbeta,nmeas;
double kappa,lambda;
double beta[nbetamax];

// Observables per measurment
double u[nmeasmax],phi2[nmeasmax],llink[nmeasmax];

// Linear contributions to the second derivative
double clin[nmeasmax],sphi2lin[nmeasmax];

// Average value of observables
double uaver[nbetamax],uerr[nbetamax]; // derivative w.r.t. beta
double caver[nbetamax],ccerr[nbetamax]; // 2nd derivative w.r.t. beta
double phi2aver[nbetamax],phi2err[nbetamax];  // derivative w.r.t kappa
double sphi2aver[nbetamax],sphi2err[nbetamax];  // 2nd. derivative w.r.t kappa
double linkaver[nbetamax],linkerr[nbetamax];  // link occupation number
double slinkaver[nbetamax],slinkerr[nbetamax];  // susceptibility of link occup. number

fstream file;

void processdata( int &argc, char *argv[] );
void finalize(int it);

char texthelp[]="Usage: exec -f [FILE]\n"
    "Analysis program of the U(1) gauge-Higgs model in the dual rep.\n"
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
  for ( int ibeta = 0; ibeta<=nbeta ; ibeta++ )
  {
    file << beta[ibeta]       << " "
         << uaver[ibeta]      << " " << uerr[ibeta]     << " "
         << caver[ibeta]      << " " << ccerr[ibeta]    << " "
         << phi2aver[ibeta]   << " " << phi2err[ibeta]  << " "
         << sphi2aver[ibeta]  << " " << sphi2err[ibeta] << " "
         << linkaver[ibeta]   << " " << linkerr[ibeta]  << " "
         << slinkaver[ibeta]  << " " << slinkerr[ibeta] << " "
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
  double efblock, lfblock, mfblock;

  double msum = 0.0;
  double esum = 0.0;
  double lsum = 0.0;
  for(im=0 ; im<nmeas ; im++)
  {
    esum = esum + u[im];
    msum = msum + phi2[im];
    lsum = lsum + llink[im];
  }
  uaver[it] = esum/(double)nmeas;
  phi2aver[it] = msum/(double)nmeas;
  linkaver[it] = lsum/(double)nmeas;

  double efsum = 0.0;
  double mfsum = 0.0;
  double lfsum = 0.0;
  for(im=0 ; im<nmeas ; im++){
    efsum += pow( uaver[it] - u[im] , 2 );
    mfsum += pow( phi2aver[it] - phi2[im] , 2 );
    lfsum += pow( linkaver[it] - llink[im] , 2 );
  }
  uerr[it] = sqrt(efsum)/(double)nmeas;
  phi2err[it] = sqrt(mfsum)/(double)nmeas;
  linkerr[it] = sqrt(lfsum)/(double)nmeas;

  for(im=0 ; im<nmeas ; im++){
    efsum += clin[im];
    mfsum += sphi2lin[im];
  }
  caver[it] = efsum/(double)nmeas;
  sphi2aver[it] = mfsum/(double)nmeas;
  slinkaver[it] = lfsum/(double)nmeas;

  double cfsum = 0.0;
  double sfsum = 0.0;
  double slsum = 0.0;
  for(im=0 ; im<nmeas ; im++){  
    efblock  = ( efsum - pow( uaver[it] - u[im], 2 ) - clin[im] )/double(nmeas-1);
    mfblock  = ( mfsum - pow( phi2aver[it] - phi2[im], 2 ) - sphi2lin[im] )/double(nmeas-1);
    lfblock  = ( lfsum - pow( linkaver[it] - llink[im], 2 ) )/double(nmeas-1);

    cfsum  += pow( caver[it] - efblock , 2 );
    sfsum  += pow( sphi2aver[it] - mfblock , 2 );
    slsum  += pow( slinkaver[it] - lfblock , 2 );
  }  
  ccerr[it] = sqrt(cfsum);
  sphi2err[it] = sqrt(sfsum);
  slinkerr[it] = sqrt(slsum);
   
  phi2aver[it]  /= (double)nsite;
  phi2err[it]   /= (double)nsite;
  sphi2aver[it] /= (double)nsite;
  sphi2err[it]  /= (double)nsite;

  uaver[it]  /= (double)(6*nsite);
  uerr[it]   /= (double)(6*nsite);
  caver[it]  /= (double)(6*nsite);
  ccerr[it]  /= (double)(6*nsite);

  linkaver[it]  /= (double)(4*nsite);
  linkerr[it]   /= (double)(4*nsite);
  slinkaver[it] /= (double)(4*nsite);
  slinkerr[it]  /= (double)(4*nsite);
}
   
//!------------------------------------------------------------------------
void processdata( int &argc, char *argv[] )
{  
  int    nequi,nskip,iseed;
  double beta0,dbeta;
  
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
  file.read((char*)&leng_t,sizeof(int));
  file.read((char*)&kappa,sizeof(double));
  file.read((char*)&lambda,sizeof(double));
  file.read((char*)&beta0,sizeof(double));
  file.read((char*)&dbeta,sizeof(double));
  file.read((char*)&nbeta,sizeof(int));
  file.read((char*)&nequi,sizeof(int));
  file.read((char*)&nmeas,sizeof(int));
  file.read((char*)&nskip,sizeof(int));
  file.read((char*)&iseed,sizeof(int));
     
  printf(" leng   = %d\n", leng);
  printf(" leng_t = %d\n", leng_t);
  printf(" kappa  = %f\n", kappa);
  printf(" lambda = %f\n", lambda);
  printf(" beta0  = %f\n", beta0);
  printf(" dbeta  = %f\n", dbeta);
  printf(" nbeta  = %d\n", nbeta);
  printf(" nequi  = %d\n", nequi);
  printf(" nmeas  = %d\n", nmeas);
  printf(" nskip  = %d\n", nskip);
  printf(" iseed  = %d\n", iseed);//*/
  
  nsite  = leng*leng*leng*leng_t;
  calculate_Pn( lambda, kappa );

  int ibeta, nblock[400], nbytes;
  for(ibeta=0 ; ibeta<=nbeta ; ibeta++ )
  {
    beta[ibeta] = beta0 + ibeta*dbeta;
    calculate_Ibeta( beta[ibeta] );

    int imeas;
    for ( imeas = 0; imeas<nmeas; imeas++ )
    {
       file.read( (char*)&nbytes, sizeof(int) ); //cout << nbytes << endl;
       file.read( (char*)nblock, nbytes*sizeof(int) );
       int inp=0, inf=0, inl=0, is_plaq=1, is_flux=0, next=1;
       u[imeas] = 0.0;
       clin[imeas] = 0.0;
       phi2[imeas] = 0.0;
       sphi2lin[imeas] = 0.0;
       llink[imeas] = 0.0;
       int total_np=0 , total_nf=0, total_nl=0;
       do{
         if(is_plaq)
         {
           int np = nblock[inp]; //cout << "np["<< inp << "] = " << np << endl;
           if ( np != -1 )
           {
             u[imeas] += np * dIn[inp];
             clin[imeas] += np * ( d2In[inp] - pow(dIn[inp], 2) );
             total_np += np;
             inp++;
           }
           else
           {
             is_plaq=0;
             is_flux=1;
           }
         }
         else if (is_flux)
         {
           int nf = nblock[inp+1+inf]; //cout << "nf["<< inf << "] = " << nf << endl;
           if ( nf != -2 )
           {
             phi2[imeas] += nf * ( Pn[inf+2]/Pn[inf] );
             sphi2lin[imeas] += nf * ( Pn[inf+4]/Pn[inf] - pow( Pn[inf+2]/Pn[inf], 2 ) );
             total_nf += nf;
             inf++;
           }
           else
           {
             is_flux=0;
           }
         }
         else
         {
           int nl = nblock[inp+inf+2+inl]; //cout << "nl["<< inl << "] = " << nl << endl;
           if ( nl != -3 )
           {
             llink[imeas] += nl*inl;
             total_nl += nl;
             inl++;
           }
           else
           {
             next= 0;
           }
         }
         //cout << "imeas= " << imeas << " - nbytes= " << nbytes << " " << endl; 
       }while(next);
       if (nbytes!=(inp+inf+inl+3))
       {
         cout << "ERROR reading data" << endl;// exit;
       }
       if ( total_np != 6*nsite )
       {
         cout << "ERROR: np != 6*nsite " << total_np <<"/"<< 6*nsite << endl; //exit;
       }
       if ( total_nf != nsite )
       {
         cout << "ERROR: nf != nsite " << total_nf << "/" << nsite << endl; //exit;
       }
       if ( total_nl != 4*nsite )
       {
         cout << "ERROR: nl != 4*nsite " << total_nl << "/" << 4*nsite <<  endl;
       }
       //else
       //{
       //  cout << imeas << ": " << nbytes << "B - np= " << inp << " - nf= " << inf << endl; 
       //}
    }
    printf("%d reading beta = %f\n", ibeta, beta[ibeta]);
    
    finalize(ibeta);
    
  }
  file.close();
}
