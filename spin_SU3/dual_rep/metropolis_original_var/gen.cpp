/************************************************************************
* Program to generate configurations for different values of mu 
* in the dual representation of the SU(3) spin model (before changing variables).  
* See ref.: arXiv:1204.6074
*
* Input file: ./bin/metro_su3_mu.start
* ./bin/metro_su3_mu.start contains the input parameters
* 
* To execute: ./bin/gen$(SIZE)_mu.x
*
* By: Ydalia Delgado (ydelgado83@gmail.com)
**************************************************************************/
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
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
    eta = kappa*exp(mu);
    etabar = kappa*exp(-mu);
    init_vtau();
    init_veta();
    init_vetabar();
    printf("%d - tau = %8.6f - kappa = %8.6f - mu = %8.6f\n",ipar,tau,kappa,mu);

    // Equilibration steps
    nsweeps(nequi);
    for (int imeas=0; imeas<nmeas; imeas++){
      // discarded MC stes
      nsweeps(nskip);
      // measurement
      measure(imeas);
    }
    file.write((char*)nblock,3*nmeas*sizeof(int));
  }
  file.close();
  rm_arrays();

  cout << "Done!!" << endl;
}

//_________________________________________________________________________
void measure(int im)
{
  im *= 3;

  nblock[im]   = ndim;
  nblock[im+1] = nmon[0];
  nblock[im+2] = nmon[1];
}


