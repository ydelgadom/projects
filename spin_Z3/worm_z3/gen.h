#ifndef _GEN_H
#define _GEN_H

//-----------------------------
// VARIABLES
//-----------------------------

char outfile[150];
fstream file;

// Variables to measure total occupation numbers
int *nblock;
int nmon[3],ndim[3];

// MC parameters
int nequi; // equilibration steps
int nmeas; // number of measurements
int nskip; // discarded steps between measurements
int iseed; // seed for Random Number Generator
int npar;   // number of parameters (tau, kappa or mu) (observables computed as a function of par)
double par0; // initial parameter
double dpar; // step size

// Model parameters
double kappa; // inverse mass
double mu; // chemical potential
double tau; // temperature

// Flux variables
int mon[nsite]; // site occupation number
int dim[nsite][3]; // link occupation number

// Weights
double monoweight[3]; // monomer weights (it depends on mu and kappa)
double bb; // weight that depends of the temperature

// Triality sum
int triadd[3][3];

//-------------------------------
// SUBROUTINES
//-------------------------------

void measure(int im);
bool check();

#endif
