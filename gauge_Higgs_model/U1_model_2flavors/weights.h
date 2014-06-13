#ifndef _WEIGHTS_H
#define _WEIGHTS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "gsl/gsl_integration.h"
#include <boost/math/special_functions/bessel.hpp>

// Variables
double In[LENGTH], Pn[LENGTH], fac[LENGTH];

#ifdef ANALYSIS
double dIn[LENGTH], d2In[LENGTH];
#endif

//______________________________________________________________________________
double fact( int q )
{
  double f = 1.0;
  for ( int j = 2; j<=q; j++) f *= j;

  return f;
}


//_________________________________________________________________________
double pn( double x, void * params )
{
  double * prms = (double *)params;
  return pow( (1.-x)/x, (int)(prms[0]+1) ) * exp( -prms[1]*pow((1.-x)/x, 2) - prms[2]*pow((1.-x)/x, 4) ) / pow( x, 2 );
}


//_________________________________________________________________________
void calculate_Pn( double ll, double kk )
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result,error;
  double params[3];
     
  gsl_function F;
  F.function = &pn;
  F.params = params;

  int i;
  params[1] = kk;
  params[2] = ll;

  gsl_set_error_handler_off();
  for ( i=0; i<LENGTH; i++)
  {
    params[0] = (double)i;
    if ( gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &result, &error) !=0 )
    {
      cout << "*** W-INTEGRATION FAILED: i=" << i << " k=" << kk << " l=" << ll << " ***" << endl;
#ifndef ANALYSIS
      MPI_Abort( MPI_COMM_WORLD, 0 );
#endif
    }
    gsl_set_error_handler(NULL);
    Pn[i] = result;
#ifndef ANALYSIS
    if ( Pn[i]>0 ){
      Pn[i] = log( Pn[i] );}
#endif
  }

  gsl_integration_workspace_free(w);
}


//_________________________________________________________________________
void calculate_Ibeta( double bb )
{
  int i;
  for ( i=0; i<LENGTH; i++)
  {
    In[i] = boost::math::cyl_bessel_i(i, bb); // gsl_sf_bessel_In( i, bb );

#ifdef ANALYSIS
    double Inp1 = boost::math::cyl_bessel_i(i+1, bb); //gsl_sf_bessel_In ( i+1, bb );
    dIn[i]  = ((double)i)/bb + Inp1/In[i];
    d2In[i] = ((double)(i*(i-1))) / pow(bb,2) + 
              ( (((double)(2*i+1))/bb) * Inp1 + gsl_sf_bessel_In(i+2, bb) ) / In[i];
#else
    if (In[i]>0) 
      In[i] = log( In[i] );
#endif
  }
}

#endif
