#ifndef _WEIGHTS_H
#define _WEIGHTS_H

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

// Variables
double In[LENGTH], Pn[LENGTH], fac[LENGTH];

#ifdef ANALYSIS
double dIn[LENGTH], d2In[LENGTH];
#endif

// Subroutines
static inline double fact( int q )
{
  double f = 1.0;
  int j;
  for (j = 2; j<=q; j++) f *= j;

  return f;
}

//_________________________________________________________________________
static double pn( double x, void * params )
{
  double * prms = (double *)params;
  double x2 = x*x;
  double result= pow( x, (int)(prms[0]+1) ) * exp( -prms[1]*x2 - prms[2]*x2*x2 );

  return result;
}

//_________________________________________________________________________
void calculate_Pn( double ll, double kk )
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
       
  double result, error;
  double params[3];
     
  gsl_function F;
  F.function = &pn;
  F.params = params;

  int i;
  params[1] = kk;
  params[2] = ll;
  for ( i=0; i<LENGTH; i++)
  {
    params[0] = (double)i;
    gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, w, &result, &error);
    Pn[i] = result;
#ifndef ANALYSIS
    if ( Pn[i]>0 )
      Pn[i] = log( Pn[i] );
#endif
  }
  gsl_integration_workspace_free (w);
}

//_________________________________________________________________________
void calculate_Ibeta( double bb )
{
  for (int i=0; i<LENGTH; i++)
  {
    In[i] = gsl_sf_bessel_In ( i, bb );
    //cout << "In " << i << " " << In[i] << endl;

#ifdef ANALYSIS
	  double Inp1 = gsl_sf_bessel_In ( i+1, bb );
    dIn[i]  = ((double)i)/bb + Inp1/In[i];
    d2In[i] = ((double)(i*(i-1))) / pow(bb,2) + 
              ( (((double)(2*i+1))/bb) * Inp1 + gsl_sf_bessel_In(i+2, bb) ) 
							/ In[i];
#else
    if (In[i]>0) 
      In[i] = log( In[i] );
#endif
  }
}

#endif
