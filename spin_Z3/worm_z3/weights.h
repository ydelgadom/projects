#ifndef _WEIGHTS_H
#define _WEIGHTS_H

//_________________________________________________________________________
inline void calcdimerweights()
{
	double expt = exp(tau);
  bb = (expt*expt - 1./expt)/(expt*expt + 2./expt);
}

//_________________________________________________________________________
inline void calcmonoweights()
{

  int l;
  const float pi2 = 2.0*3.141592653589793;
    
  for (l = -1; l<=1; ++l){  
    monoweight[l+1] = ( exp(2.0*kappa*cosh(mu)) + 2.0*exp(-kappa*cosh(mu))*
                      cos(sqrt(3.0)*kappa*sinh(mu) - l*pi2/3.0) )/3.0;

    cout << "M["<< l <<"] = " << monoweight[l+1] << endl; 
  }
}

#endif
