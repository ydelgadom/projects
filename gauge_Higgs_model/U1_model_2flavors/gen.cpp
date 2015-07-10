/************************************************************************
* Program to generate configurations for different values of beta 
* in the dual representation of the U(1) gauge-Higgs model (2 flavors).
* See references in ../papers
* 
* To execute: 
*     ./bin/gen$(SIZE).x -[parameters] (list of input parameters in init.h)
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
#include <sstream>
#include <getopt.h>

using namespace std;

extern "C"
{
#include "ranlxd.h"
}

#define iseed 12245

#include "lattice.h"
#include "gen.h"
#include "weights.h"
#include "init.h"
#include "multigrid.h"
#include "worm.h"
#include "sweeps.h"

//______________________________________________________________________________

int main(int argc, char *argv[])
{
  //  Initialize MPI.
  MPI::Init ( argc, argv );

  //  Get the number of processes.
  nproc = MPI::COMM_WORLD.Get_size ( );

  //  Determine this process's rank.
  id = MPI::COMM_WORLD.Get_rank ( );

  //  Process input parameters and initilizes arrays,lattice, etc...
  init( argc, argv, iseed + 100*id );

  //  Each core generates configurations for 1 set of parameters
  beta = beta0 + id*dbeta;
  mu = mu0 + id*dmu;
  kappa = kappa0 + id*dkappa;

  //  Compute weights
  calculate_Ibeta( beta );
  calculate_Pn( lambda, kappa );

  int gslerr;
  if ( sub_leng != leng && sub_leng>0 )
  {
    cout << "Multigrid equilibration " << endl;
    // 90% of equilibration steps performed in a smaller lattice.
    gslerr = multigrid_equi( );

    // 10% of equilibration steps in the full lattice
    gslerr = nsweeps( nequi/10, leng, nsite, vneib );
  }
  else
  {
    cout << "Normal equilibration " << endl;

    // all equilibration steps in the full lattice
    gslerr = nsweeps( nequi, leng, nsite, vneib );
  }
  for (int imeas=0 ; imeas<nmeas ; imeas++ )
  { 
    // do nskip steps between 2 measurements
    gslerr = nsweeps( nskip, leng, nsite, vneib );
    int nbytes = measure( gslerr );

    // write to disk after each measurement.
    file.write( (char*)nblock, nbytes );
  }
  if ( gslerr != -1 ) write_conf( );
#ifdef CHECK
  check( nsite, vneib );
#endif
  file.close( );

  printf( "%d - (b, k, m) = (%5.4f, %4.3f, %4.3f)\n\n", id, beta, kappa, mu );

  //  Terminate MPI.
  MPI::Finalize( );

  return 0;
}

//_________________________________________________________________________
int nsweeps( int ns, int vleng, int lnsite, int neib[][8] )
{
  for( int is=1; is<=ns; ++is )
  {
#ifdef SWA
    int inw = 0;
    while ( inw < nsite )
    {
      ++inw;
      if ( surface_worms( inw%2, lnsite, neib ) < 0 ) return -1;
    };
#endif

    if ( sweep_lprime( lnsite, neib ) < 0 ) return -1;

#ifdef CUBES
    sweep_cube( );
#endif

#ifdef WL
    if ( sweep_winding_loop( vleng, vneib ) < 0 ) return -1;
#endif
  }
  return 0;
}

//_____________________________________________________________________________
int measure( int gslerr )
{
  if ( gslerr != -1 )
  {
    int i_pmax=0, i_fmax[2]={0,0};
    int i;
    for ( i=0; i<LENGF; ++i )
    {
      if (nplaq[i]) i_pmax = i;
      if (nflux[0][i]) i_fmax[0] = i;
      if (nflux[1][i]) i_fmax[1] = i;
    }
    ++i_pmax;
    ++i_fmax[0];
    ++i_fmax[1];

    i = 3;
    memcpy( &nblock[i], nplaq, i_pmax*sizeof(int) );
    i += i_pmax;
    nblock[i] = -1; ++i;

    memcpy( &nblock[i], &nflux[0][0], i_fmax[0]*sizeof(int) );
    i += i_fmax[0]; 
    nblock[i] = -2; ++i;

    memcpy( &nblock[i], &nflux[1][0], i_fmax[1]*sizeof(int) );
    i += i_fmax[1];
    nblock[i] = -3; ++i;

    nblock[0] = i - 1;
    nblock[1] = nlink[0];
    nblock[2] = nlink[1];
  
    return i*sizeof(int);
  }
  else
  {
    nblock[0] = 8;
    nblock[1] = 0;
    nblock[2] = 0;
    nblock[3] = 6*nsite;
    nblock[4] = -1;
    nblock[5] = nsite;
    nblock[6] = -2;
    nblock[7] = nsite;
    nblock[8] = -3;

    return 9*sizeof(int);
  }
}

//_________________________________________________________________________
void write_conf( )
{
  string ss;
  std::ostringstream oss(ss);
  oss << outfile << ".conf";
  fconf.open( oss.str( ).c_str( ),ios::trunc | ios::out | ios::binary );
  fconf.write( (char*)&vlink[0][0], 8*nsite*sizeof(int) );
  fconf.write( (char*)&vllink[0][0], 8*nsite*sizeof(int) );
  fconf.write( (char*)&vplaq[0][0], 6*nsite*sizeof(int) );
  fconf.write( (char*)&vflux[0][0], 2*nsite*sizeof(int) );
  fconf.close( );
}

#ifdef CHECK
//_______________________________________________________________________
void check( int lnsite, int neib[][8] )
{
  int error = 0;
  int serror = 0;
  int is;
  for ( is=0; is<lnsite; is++ )
  {
    int val = vlink[is][0]        + vlink[is][1] 
          + vlink[is][2]          + vlink[is][3]
          - vlink[neib[is][4]][0] - vlink[neib[is][5]][1] 
          - vlink[neib[is][6]][2] - vlink[neib[is][7]][3];
    if (val!=0) serror++;//cout << "ERROR: site " << is << " ; sume = " << val << endl;

    val = vlink[is][4]        + vlink[is][5] 
          + vlink[is][6]      + vlink[is][7]
          - vlink[neib[is][4]][4] - vlink[neib[is][5]][5] 
          - vlink[neib[is][6]][6] - vlink[neib[is][7]][7];
    if (val!=0) serror++;//cout << "ERROR: site " << is << " ; sume = " << val << endl;

    val = vlink[is][0] - vlink[is][4]  
          + vplaq[is][0+1-1] - vplaq[neib[is][5]][0+1-1] 
          + vplaq[is][0+2-1] - vplaq[neib[is][6]][0+2-1] 
          + vplaq[is][0+3  ] - vplaq[neib[is][7]][0+3  ];
    if (val!=0) error++;//cout << "ERROR: link 0 " << is << " ; sume = " << val << endl;

    val = vlink[is][1] - vlink[is][5]
          - vplaq[is][1+0-1] + vplaq[neib[is][4]][1+0-1] 
          + vplaq[is][1+2-1] - vplaq[neib[is][6]][1+2-1] 
          + vplaq[is][1+3  ] - vplaq[neib[is][7]][1+3  ];
    if (val!=0) error++;//cout << "ERROR: link 1 " << is << " ; sume = " << val << endl;

    val = vlink[is][2] - vlink[is][6]
          - vplaq[is][2+0-1] + vplaq[neib[is][4]][2+0-1] 
          - vplaq[is][2+1-1] + vplaq[neib[is][5]][2+1-1] 
          + vplaq[is][2+3  ] - vplaq[neib[is][7]][2+3  ];
    if (val!=0) error++;//cout << "ERROR: link 2 " << is << " ; sume = " << val << endl;

    val = vlink[is][3] - vlink[is][7]
          - vplaq[is][3+0] + vplaq[neib[is][4]][3+0] 
          - vplaq[is][3+1] + vplaq[neib[is][5]][3+1] 
          - vplaq[is][3+2] + vplaq[neib[is][6]][3+2];
    if (val!=0) error++;//cout << "ERROR: link 3 " << is << " ; sume = " << val << endl;
  }
  if( serror>0 ) cout << "ERROR SITE " << serror << endl;
  if(  error>0 ) cout << "ERROR LINK " <<  error << endl;
}
#endif
