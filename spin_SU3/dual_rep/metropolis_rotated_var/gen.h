#ifndef _GEN_H
#define _GEN_H

//-----------------------------
// VARIABLES
//-----------------------------

#define MAXFLUX 1

char outfile[150];
fstream file;

// Arrays for measurements
int *nblock;
int ndim,nmon[3];

// MC parameters
int nequi; // number of equilibration steps
int nmeas; // number of measurements
int nskip; // number of discarded steps between measurements
int iseed; // seed for random number generator

// Model parameters
double kappa; // inverse mass
double mu,emu; // chemical potentail and fugacity term
double tau; // temperature
double dpar,par0;
int npar;

// Flux variables per site/link
int mon[nsite][2]; // site occupation number
int dim[nsite][6]; // link occupation number
int mnx[nsite][2]; // total flux

// Weights to compute the Metropolis probability
double Tmn[MAXFLUX][MAXFLUX],facn[MAXFLUX];
double vtau[15],vkappa[15],vemu[15];
int    vnu[3] = {1,2,0};

//-------------------------------
// SUBROUTINES
//-------------------------------

void measure(int im);

#endif
