/************************************************************************
* Program to generate configurations for different values of beta 
* in the dual representation of the U(1) gauge-Higgs model (only 1 flavor).
* See references in ../papers
*
* Input file: ./bin/worm_beta.start
* ./bin/worm_beta.start contains the input parameters
* 
* To execute: ./bin/gen$(SIZE).x
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
using namespace std;

extern "C"
{
#include "ranlxd.h"
}
#include "lattice.h"
#include "gen.h"
#include "weights.h"
#include "sweeps.h"
#include "worm.h"
#include "init.h"

// MAIN PROGRAM __________________________________________________

int main( )
{
  printf("\nProgram gen.cpp\n\nCompiled with leng = %d\n",leng);

  init( );
  calculate_Pn( lambda, kappa );

  for (int ibeta=0; ibeta<=nbeta; ibeta++)
  {
    beta = beta0 + ibeta*dbeta;
    printf("%d - beta = %8.6f\n", ibeta, beta);

    calculate_Ibeta( beta );

    nworms( nequi );
    for (int imeas=0; imeas<nmeas; imeas++)
    {
      nworms( nskip ); 
      int nbytes = measure( );
      file.write( (char*)nblock, nbytes );
    }
    //check( );
  }
  file.close();
  cout << "Done." << endl;

  return 0;
}

//_________________________________________________________________________
//_________________________________________________________________________
int measure( )
{
  int i_pmax=0, i_fmax=0, i_lmax=0;
  int i;
  for ( i=0; i<LENGTH; i++ )
    if (nplaq[i]) i_pmax = i;

  for ( i=0; i<LENGTH; i++ )
    if (nlink[i]) i_lmax = i;
  
  for ( i=0; i<LENGTH_FLUX; i++ )
    if (nflux[i]) i_fmax = i;

  nblock[0] = i_pmax + i_fmax + i_lmax + 6; 

  for ( i=0; i<=i_pmax; i++)
    nblock[i+1] = nplaq[i];
  nblock[i_pmax+2] = -1;

  for ( i=0; i<=i_fmax; i++)
    nblock[i_pmax+3+i] = nflux[i];
  nblock[i_pmax+i_fmax+4] = -2;

  for ( i=0; i<=i_lmax; i++)
    nblock[i_pmax+i_fmax+5+i] = nlink[i];
  nblock[i_pmax+i_fmax+i_lmax+6] = -3;

  return (i_pmax + i_fmax + i_lmax + 7)*sizeof(int);
}

//_______________________________________________________________________
void check()
{
  int error = 0;
  int serror = 0;
  int is;
  for ( is=0; is<nsite; is++ )
  {
    int val = vlink[is][0]        + vlink[is][1] 
          + vlink[is][2]          + vlink[is][3]
          - vlink[neib[is][4]][0] - vlink[neib[is][5]][1] 
          - vlink[neib[is][6]][2] - vlink[neib[is][7]][3];

    if (val!=0) serror++;//cout << "ERROR: site " << is << " ; sume = " << val << endl;

    int il=0;
    val = vlink[is][il]  
          + vplaq[is][0+1-1] - vplaq[neib[is][5]][0+1-1] 
          + vplaq[is][0+2-1] - vplaq[neib[is][6]][0+2-1] 
          + vplaq[is][0+3  ] - vplaq[neib[is][7]][0+3  ];
    if (val!=0) error++;//cout << "ERROR: link 0 " << is << " ; sume = " << val << endl;

    il=1;
    val = vlink[is][il]
          - vplaq[is][1+0-1] + vplaq[neib[is][4]][1+0-1] 
          + vplaq[is][1+2-1] - vplaq[neib[is][6]][1+2-1] 
          + vplaq[is][1+3  ] - vplaq[neib[is][7]][1+3  ];
    if (val!=0) error++;//cout << "ERROR: link 1 " << is << " ; sume = " << val << endl;

    il=2;
    val = vlink[is][il]
          - vplaq[is][2+0-1] + vplaq[neib[is][4]][2+0-1] 
          - vplaq[is][2+1-1] + vplaq[neib[is][5]][2+1-1] 
          + vplaq[is][2+3  ] - vplaq[neib[is][7]][2+3  ];
    if (val!=0) error++;//cout << "ERROR: link 2 " << is << " ; sume = " << val << endl;

    il=3;
    val = vlink[is][il]
          - vplaq[is][3+0] + vplaq[neib[is][4]][3+0] 
          - vplaq[is][3+1] + vplaq[neib[is][5]][3+1] 
          - vplaq[is][3+2] + vplaq[neib[is][6]][3+2];
    if (val!=0) error++;//cout << "ERROR: link 3 " << is << " ; sume = " << val << endl;
  }
  if( serror>0 ) cout << "ERROR SITE " << serror << endl;
  if(  error>0 ) cout << "ERROR LINK " <<  error << endl;
}
