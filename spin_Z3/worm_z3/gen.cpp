/************************************************************************
* Program to generate configurations for different values of $PAR (=TAU,MU or KAPPA)
* in the dual representation of the Z(3) spin model.
* The program can also use open worms (flag -DOPEN during compilation).
* See references in ../papers
*
* Input file: ./bin/worm_$PAR.start
* ./bin/worm_$PAR.start contains the input parameters
*
* To compile
*   closed worms: make SIZE=4 PAR=TAU gen
*   open worms  : make SIZE=4 PAR=TAU genopen
* 
* To execute: ./bin/gen$(SIZE)_$PAR.x
*         or: ./bin/gen$(SIZE)_open_$PAR.x
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <ctime>
using namespace std;

extern "C"
{
#include "ranlxd.h"
}

#include "lattice.h"
#include "gen.h"
#include "init.h"
#include "weights.h"
#include "worm.h"

// MAIN PROGRAM __________________________________________________

int main()
{
  printf("\nProgram gen.cpp\n\nCompiled with leng = %d\n",leng);
#ifdef OPEN
  printf("Compiled for open worms\n");
#else
  printf("Compiled for closed worms\n");
#endif

  read_params();
  rlxd_init(1,iseed);
  init_fields(); 
  mk_arrays();  
    
  for ( int ipar=0; ipar<=npar; ++ipar )
  {
#ifdef KAPPA
    kappa = par0 + dpar*ipar;
#endif
#ifdef MU
    mu = par0 + dpar*ipar;
#endif
#ifdef TAU
    tau = par0 + dpar*ipar;
#endif
    printf("ipar = %d - tau= %8.6f - mu= %8.6f - kappa = %8.6f\n",ipar,tau,mu,kappa);

    // Calculate monomer and dimer weights 
    calcdimerweights();
    calcmonoweights();

    // Thermalization steps
    nworms(nequi);

    // Measurements
    for (int imeas = 0; imeas<nmeas; imeas++)
    {
      // discarded steps
      nworms(nskip);

      // measurement
      measure(imeas);
    }
    file.write((char*)nblock,6*nmeas*sizeof(int));
  }
  file.close();
  rm_arrays();

  cout << "Done!!!" << endl;

}

//_________________________________________________________________________
void measure(int im)
{
  nblock[im*6+0] = ndim[0];
  nblock[im*6+1] = ndim[1];
  nblock[im*6+2] = ndim[2];

  nblock[im*6+3] = nmon[0];
  nblock[im*6+4] = nmon[1];
  nblock[im*6+5] = nmon[2];
}

//_________________________________________________________________________
bool check()
{
  int is,isn,inu,dimsum;

  for(is=0;is<nsite;++is){
    dimsum = 0;
    for (inu=0;inu<3;inu++){
      dimsum += (dim[is][inu]-1);
      isn = neib[is][inu+3];
      dimsum -= (dim[isn][inu]-1);
    }

    dimsum += (mon[is]-1);
    if (dimsum%3 != 0 ){ cout << "WRONG!!!!" << endl; return false; }

  }
  return true;
}

