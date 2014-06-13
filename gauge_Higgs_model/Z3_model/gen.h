#ifndef _GEN_H
#define _GEN_H

/******************************
* Global variables
******************************/

fstream file;
char outfile[150];

/* average occupation numbers that will 
*  be written on the output file */
int *nblock;

/* varialbes with occupation number */
int vplaq[nsite][6];
int vlink[nsite][4];

/* average occupation numbers */
int nlink[4][3],nplaq[3];

/* parameters of the MC program */
int nskip;  // steps for decorrelation
int nequi;  // thermalization steps
int nmeas;  // number of measurements
int iseed;  // seed for the random number generator
double beta;   // coupling constant of the gauge action 
double gama;   // parameter of the Higgs field
double mu;     // current value of the chemical potential
double par0,dpar;
int npar;

/********************************
* Subroutines
********************************/

void measure( int im );
void check();   // checks if worm update is correct

#endif
