#ifndef _GEN_H
#define _GEN_H

#define SWA // generate configs. using SWA
#define WL // use winding loop update
//#define CUBES // use sweeps of cube update
//#define CHECK // check if updates are correct

#define LENGTH 50
#define LENGF 50

#include "mpi.h"

//--------------------------------------------------
// Global Variables
//--------------------------------------------------

// Files:
fstream file, fconf;
char outfile[128]; // output with measurements.
char conffile[128]; // file where initial config. is stored.

// MC parameters
int nskip; // discarded steps between measurements
int nequi; // equilibration steps 
int nmeas; // measurements steps

double lambda; 
double beta0, dbeta, beta;
double kappa0, dkappa, kappa;
double mu0, dmu, mu;

// Flag to read a initial config.:
// rconf = false ==> cold start
// rconf = true ==> read config. and start with it. 
bool rconf;

// Constrained occupation number: 
// vlink[][0...3] = first flavor, vlink[][4...7] = second flavor
int vlink[nsite][8];

// Unconstrained occupation number: 
// vllink[][0...3] = first flavor, vllink[][4...7] = second flavor
int vllink[nsite][8];

// Plaquette occupation number
int vplaq[nsite][6];

// Total sum of links and plaquette occupation numbers per site.
// vflux[][0] = first flavor, vflux[][1] = second flavor
int vflux[nsite][2];

// MPI variables: id=local id and 
// nproc = total number of cores in the job
int id, nproc;

// Variables for measurements
int nblock[3*LENGF+3];
int nplaq[LENGF], nflux[2][LENGF];
int nlink[2];

//-------------------------------
// Subroutines
//-------------------------------

int nsweeps( int ns, int vleng, int lnsite, int neib[][8] );
int measure( int gslerr );
void write_conf( );
#ifdef CHECK
void check( int lnsite, int neib[][8] );
#endif

#endif 
