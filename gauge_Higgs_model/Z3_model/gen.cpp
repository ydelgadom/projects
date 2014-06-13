/************************************************************************
* Program to generate configurations for different values of $PAR (BETA or MU)
*	in the dual representation of the Z(3) gauge-Higgs model (only 1 flavor).
*	See references in ../papers
*
*	Input file: ./bin/worm_$PAR.start
* ./bin/worm_$PAR.start contains the input parameters
*	
*	To execute: ./bin/gen$(SIZE)_$PAR.x
*
*	By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
using namespace std;

extern "C"
{
/* random number generator */
#include "ranlxd.h"
}
#include "lattice.h"
#include "gen.h"
#include "weights.h"
#include "worm.h"
#include "init.h"

//____________________________________________________________________________
/* MAIN PROGRAM */

int main( )
{
  printf("\nProgram gen.cpp\n\nCompiled with lengths: %d, %d\n",leng, leng_t);

  init( );
  for (int ipar=0; ipar<=npar; ipar++)
  {
#ifdef MU
    mu = par0 + ipar*dpar;
#else
		beta = par0 + ipar*dpar;
#endif
    printf( "%d - beta = %8.6f - mu = %8.6f\n", ipar, beta, mu );

		/* compute link and site weights*/
	  calculate_bbt_weight( );
    calculate_link4weights();

    /* generates nequi*nsite worms for equilibration */
    nworms( nequi );
    for (int imeas=0; imeas<nmeas; imeas++)
    {
      /* between 2 consecutive measurements 
         nskip*nsite worms are generated for 
         decorrelation */
      nworms( nskip );
      measure( imeas );
    }
    // check();
    file.write((char*)nblock,9*nmeas*sizeof(int));
  }
  file.close( );
  /* delete nblock array */
  rm_arrays( );

  cout << "Done." << endl;
}

//_________________________________________________________________________
void measure( int im )
{
  /* function to measure the average occupation number of the elements */
  im *= 9;

  nblock[im+0] = nlink[0][0]+nlink[1][0]+nlink[2][0];
  nblock[im+1] = nlink[0][1]+nlink[1][1]+nlink[2][1];
  nblock[im+2] = nlink[0][2]+nlink[1][2]+nlink[2][2];

  nblock[im+3] = nlink[3][0];
  nblock[im+4] = nlink[3][1];
  nblock[im+5] = nlink[3][2];

  nblock[im+6] = nplaq[0];
  nblock[im+7] = nplaq[1];
  nblock[im+8] = nplaq[2];
}

//_______________________________________________________________________
void check()
{
  int is,il,val,error,serror;

  error = 0;
  serror = 0;
  for ( is=0; is<nsite; is++ )
  {
    val =   vlink[is][0]          + vlink[is][1] 
          + vlink[is][2]          + vlink[is][3]
          - vlink[neib[is][4]][0] - vlink[neib[is][5]][1] 
          - vlink[neib[is][6]][2] - vlink[neib[is][7]][3];

    if (val%3!=0) serror++;//cout << "ERROR: site " << is << " ; sume = " << val << endl;

    il=0;
    val = vlink[is][il] - 4 + 3 
          + vplaq[is][0+1-1] - vplaq[neib[is][5]][0+1-1] 
          + vplaq[is][0+2-1] - vplaq[neib[is][6]][0+2-1] 
          + vplaq[is][0+3  ] - vplaq[neib[is][7]][0+3  ];
    if (val%3!=0) error++;//cout << "ERROR: link 0 " << is << " ; sume = " << val << endl;

    il=1;
    val = vlink[is][il] - 4 + 3
          - vplaq[is][1+0-1] + vplaq[neib[is][4]][1+0-1] 
          + vplaq[is][1+2-1] - vplaq[neib[is][6]][1+2-1] 
          + vplaq[is][1+3  ] - vplaq[neib[is][7]][1+3  ];
    if (val%3!=0) error++;//cout << "ERROR: link 1 " << is << " ; sume = " << val << endl;

    il=2;
    val = vlink[is][il] - 4 + 3
          - vplaq[is][2+0-1] + vplaq[neib[is][4]][2+0-1] 
          - vplaq[is][2+1-1] + vplaq[neib[is][5]][2+1-1] 
          + vplaq[is][2+3  ] - vplaq[neib[is][7]][2+3  ];
    if (val%3!=0) error++;//cout << "ERROR: link 2 " << is << " ; sume = " << val << endl;

    il=3;
    val = vlink[is][il] - 4 + 3
          - vplaq[is][3+0] + vplaq[neib[is][4]][3+0] 
          - vplaq[is][3+1] + vplaq[neib[is][5]][3+1] 
          - vplaq[is][3+2] + vplaq[neib[is][6]][3+2];
    if (val%3!=0) error++;//cout << "ERROR: link 3 " << is << " ; sume = " << val << endl;
  }
  if( serror>0 ) cout << "ERROR SITE " << serror << endl;
  if(  error>0 ) cout << "ERROR LINK " <<  error << endl;
}

