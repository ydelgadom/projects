#ifndef _WEIGHTS_H
#define _WEIGHTS_H

// Variables
/* weights for the Metropolis step during the updates */
double bb, bbt, logbb, logbbt;
double link4weight[3];

// Subroutines
//______________________________________________________________________
void inline calculate_bb_weight( )
{
  bb  = ( exp(2.*gama) - exp(-gama) )/( exp(2.*gama) + 2.*exp(-gama) );
  logbb = log(bb);
}

//______________________________________________________________________
void inline calculate_bbt_weight( )
{
  bbt = ( exp(beta) - exp(-beta/2.) )/( exp(beta) + 2.*exp(-beta/2.) );
  logbbt = log(bbt);
}

//______________________________________________________________________
void calculate_link4weights()
{
  const double pi2over3 = 2.0*3.141592653589793/3.0;

  for (int i = -1; i<=1; i++)
  {
    link4weight[i+1] = ( exp(2.*gama*cosh(mu)) + 2.*exp(-gama*cosh(mu) ) *
                         cos(sinh(mu)*sqrt(3.)*gama - i*pi2over3) )/3.0;

    cout << "M["<< i <<"] = " << link4weight[i+1] << endl;

    link4weight[i+1] = log(link4weight[i+1]);
  }
  cout << endl;
}

#endif
