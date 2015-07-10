/******************************************************************
* Generate configurations in the conventional representation of the 
* SU(3) spin model.  See ref.: arXiv:1204.6074
*
* Input file: ./bin/metro_su3.start
* 
* To execute: ./bin/gen$(SIZE).x
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
******************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <ctime>
using namespace std;

extern "C"
{
#include "ranlxd.h"
}

fstream file;
char outfile[150];

const double pi   = 3.141592653589793;
const int    leng = SIZE;
const int    nsite = leng*leng*leng;

int    neib[nsite][6];
int    nskip,nequi,nmeas,iseed,ntau;
double dtau,tau0,kappa,eps,tau;
double theta2[nsite],theta1[nsite],mm[nsite];
double px[2][nsite];

double *nblock;

void init_lattice();
void read_params();
void init_fields();
void nsweeps(int ns);
void measure(int im);

// MAIN PROGRAM __________________________________________________

int main()
{
  int itau,imeas;
  

  clock_t t1 = clock();

  printf("\nProgram worm_su3_metropolis.c\n\nCompiled with leng = %d\n",leng);

  read_params();
  rlxd_init(1,iseed);
  init_fields();

  nblock = new double[2*nmeas];
     
  for (itau=0; itau<=ntau; itau++)
  {
    tau = tau0 + itau*dtau;
    printf("%d - tau = %8.6f\n",itau,tau);

    nsweeps(nequi);  
    for (imeas=0; imeas<nmeas; imeas++)
    {
      nsweeps(nskip);
      measure(imeas);
    }
    file.write((char*)nblock,2*nmeas*sizeof(double));
  }
  file.close();

  delete[] nblock;

  cout << "Done in " << (clock()-t1) << endl;

}

//_________________________________________________________________________
//_________________________________________________________________________

void measure(int im)
{
  double pp,ee,pneib[2];

  im *= 2;

  ee = 0.;
  pp = 0.;
  for (int ix=0; ix<nsite; ix++)
  {
    pneib[0] = px[0][neib[ix][0]] + px[0][neib[ix][1]] + px[0][neib[ix][2]];
    pneib[1] = px[1][neib[ix][0]] + px[1][neib[ix][1]] + px[1][neib[ix][2]];
    
    ee += 2.*kappa*px[0][ix] + 2.*tau*( px[0][ix]*pneib[0] + px[1][ix]*pneib[1] );
    pp += px[0][ix];
  }
  nblock[im]   = -ee;
  nblock[im+1] = pp;

}
//_________________________________________________________________________

void nsweeps(int ns)
{
  int    is,ix;
  double mmnew, theta1new, theta2new, ran[3];
  double pneib[2],rho,pxnew[2],rhoold;
  
  for(is=1; is<=ns; is++)
  {
    for(ix=0; ix<nsite; ix++)
    {
      ranlxd(ran,3);
      theta1new = theta1[ix] + 2.*eps*(ran[0]-0.5);
      theta2new = theta2[ix] + 2.*eps*(ran[1]-0.5); 

      pxnew[0] = cos(theta1new) + cos(theta2new) + cos(theta1new+theta2new);
      pxnew[1] = sin(theta1new) + sin(theta2new) - sin(theta1new+theta2new);

      mmnew = 8.*( 1. - cos(theta1new-theta2new) )*( 1. - cos(2.*theta1new+theta2new) )*( 1. - cos(theta1new+2.*theta2new) );

      pneib[0] = px[0][neib[ix][0]] + px[0][neib[ix][1]] + px[0][neib[ix][2]] + 
                 px[0][neib[ix][3]] + px[0][neib[ix][4]] + px[0][neib[ix][5]];
      pneib[1] = px[1][neib[ix][0]] + px[1][neib[ix][1]] + px[1][neib[ix][2]] + 
                 px[1][neib[ix][3]] + px[1][neib[ix][4]] + px[1][neib[ix][5]];
      rhoold = exp( 2.*kappa*(px[0][ix]) + 2.*tau*( px[0][ix]*pneib[0] + px[1][ix]*pneib[1] ) );
      rho    = exp( 2.*kappa*(pxnew[0] ) + 2.*tau*( pxnew[0]*pneib[0]  + pxnew[1]*pneib[1]  ) );
      rho    = rho*mmnew/(rhoold*mm[ix]);

      if (ran[2] < rho)
      {
        px[0][ix]  = pxnew[0];
        px[1][ix]  = pxnew[1];

        while ( theta1new > pi ){ 
          theta1new = theta1new - 2.*pi;
        };

        while ( theta1new < -pi ){ 
          theta1new = theta1new + 2.*pi;
        };

        while ( theta2new > pi ){ 
          theta2new = theta2new - 2.*pi;
        };

        while ( theta2new < -pi ){ 
          theta2new = theta2new + 2.*pi;
        };
      
        theta1[ix] = theta1new;
        theta2[ix] = theta2new;
        mm[ix]     = mmnew;
      }

    }

  }

}


//_________________________________________________________________________

void read_params()
{
  char dummy[30];

  file.open("metro_su3.start", ios::in);
  
  file >> dummy;
  file >> tau0;
  printf(" tau0    = %f\n", tau0);
  file >> dummy;
  file >> dtau;
  printf(" dtau    = %f\n", dtau);
  file >> dummy;
  file >> ntau;
  printf(" ntau    = %d\n", ntau);
  file >> dummy;
  file >> kappa;
  printf(" kappa   = %f\n", kappa);
  file >> dummy;
  file >> eps;
  printf(" eps      = %f\n", eps);
  file >> dummy;
  file >> nequi;
  printf(" nequi   = %d\n", nequi);
  file >> dummy;
  file >> nmeas;
  printf(" nmeas   = %d\n", nmeas);
  file >> dummy;
  file >> nskip;
  printf(" nskip   = %d\n", nskip);
  file >> dummy;
  file >> iseed;
  printf(" iseed   = %d\n", iseed);
  file >> dummy;
  file >> outfile;
  printf(" outfile = %s\n", outfile);
  
  file.close();
   
  file.open(outfile,ios::trunc | ios::out | ios::binary );
 
  file.write((char*)&leng,sizeof(int));
  file.write((char*)&tau0,sizeof(double));
  file.write((char*)&dtau,sizeof(double));
  file.write((char*)&ntau,sizeof(int));
  file.write((char*)&kappa,sizeof(double));
  file.write((char*)&eps,sizeof(double));
  file.write((char*)&nequi,sizeof(int));
  file.write((char*)&nmeas,sizeof(int));
  file.write((char*)&nskip,sizeof(int));
  file.write((char*)&iseed,sizeof(int));

}
//_________________________________________________________________________

void init_fields()
{
  init_lattice();

  for(int ix=0; ix<nsite; ix++ )
  {
    theta2[ix] = 0.;
    theta1[ix] = pi/2.;
    mm[ix]     = 16.;
    px[0][ix]  = 1.;
    px[1][ix]  = 0.;
  }

}
//_________________________________________________________________________

void init_lattice()
{

  int i1,i2,i3,i1p,i2p,i3p,i1m,i2m,i3m;
  int is,isp1,isp2,isp3,ism1,ism2,ism3;

  for (i1 = 0; i1<leng ; i1++ ){
    i1p = i1 + 1;
    i1m = i1 - 1;
    if (i1 == (leng-1) ) i1p = 0;
    if (i1 == 0)         i1m = leng-1;

  for (i2 = 0; i2<leng ; i2++){
    i2p = i2 + 1;
    i2m = i2 - 1;
    if (i2 == (leng-1) ) i2p = 0;
    if (i2 == 0)         i2m = leng-1;

  for (i3 = 0; i3<leng ; i3++){
    i3p = i3 + 1;
    i3m = i3 - 1;
    if (i3 == (leng-1) ) i3p = 0;
    if (i3 == 0)         i3m = leng-1;

  is = i1 + i2*leng + i3*leng*leng;  

  isp1 = i1p + i2*leng +  i3*leng*leng;   
  isp2 = i1 +  i2p*leng + i3*leng*leng; 
  isp3 = i1 +  i2*leng +  i3p*leng*leng;  
 
  ism1 = i1m + i2*leng +  i3*leng*leng;  
  ism2 = i1 +  i2m*leng + i3*leng*leng; 
  ism3 = i1 +  i2*leng +  i3m*leng*leng;

  // fill in the neighborhood array

  neib[is][0] = isp1;
  neib[is][1] = isp2;
  neib[is][2] = isp3;
 
  neib[is][3] = ism1;
  neib[is][4] = ism2;
  neib[is][5] = ism3;

  }
  }
  }

}

