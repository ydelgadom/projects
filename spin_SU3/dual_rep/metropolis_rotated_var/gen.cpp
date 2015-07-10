/************************************************************************
* Program to generate configurations for different values of $PAR (=MU, TAU) 
*	in the dual representation of the SU(3) spin model (after changing variables).  
*	See ref.: arXiv:1204.6074
*
*	Input file: ./bin/metro_su3_$PAR.start
* ./bin/metro_su3_$PAR.start contains the input parameters
*
* To compile: make SIZE=4 PAR=TAU gen
*	
*	To execute: ./bin/gen$(SIZE)_rot_$PAR.x
*
*	By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
using namespace std;

extern "C"
{
#include "ranlxd.h"
}

#include "weights.h"
#include "lattice.h"
#include "gen.h"
#include "sweeps.h"
#include "init.h"

// MAIN PROGRAM __________________________________________________

int main()
{
  printf("\nProgram gen.cpp\n\nCompiled with leng = %d\n",leng);

  read_params();
  mk_arrays();
  rlxd_init(1,iseed);
  init_fields();
     
  for (int ipar=0; ipar<=npar; ipar++)
  {
#ifdef MU
    mu = par0 + ipar*dpar;
#else
    tau = par0 + ipar*dpar;
#endif
    printf("%d - tau = %8.6f - kappa = %8.6f - mu = %8.6f\n",ipar,tau,kappa,mu);

    cout << "MAXFLUX = " << MAXFLUX << endl;

    // Compute dimer and monomer weights
    init_vtau();
    init_vemu();
    init_vkappa();

    // Thermalization steps
    try {
      nsweeps(nequi);  
      for (int imeas=0; imeas<nmeas; imeas++){
        // Discarded sweeps
        nsweeps(nskip);
        // Measurements
        measure(imeas);
      }
      file.write((char*)nblock,4*nmeas*sizeof(int));
    }
    catch (exception& e) {
       cout << e.what() << '\n';
    }
  }
  file.close();
  rm_arrays();
  
  cout << "Done!!" << endl;

}

//_________________________________________________________________________
void measure(int im)
{
  im *= 4;
  nblock[im]   = ndim;
  nblock[im+1] = nmon[0];
  nblock[im+2] = nmon[1];
  nblock[im+3] = nmon[2];
}

