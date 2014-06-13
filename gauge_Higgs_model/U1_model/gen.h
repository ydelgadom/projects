#ifndef _GEN_H
#define _GEN_H

#define LENGTH 50
#define LENGTH_FLUX 100

// Variables
fstream file;
char outfile[128];

// Flux variables
int vlink[nsite][4];
int vllink[nsite][4]; // unconstrained link
int vplaq[nsite][6];
int vflux[nsite];

// Measurements
int nblock[2*LENGTH+LENGTH_FLUX];
int nplaq[LENGTH], nflux[LENGTH_FLUX], nlink[LENGTH];

// Simulation parameters
int nskip; // discarded MC steps
int nequi; // equilibration steps
int nmeas; // number of measurements
int iseed; // seed for random number generator
int nbeta; // number of betas to be generated

// Model parameters
double lambda,kappa;
double dbeta,beta0,beta; // beta is the coupling constant

//  Subroutines
int  measure( );
void check( );

#endif
