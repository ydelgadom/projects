#ifndef _MOMENTS_H
#define _MOMENTS_H

//-------------------------------
// SUBROUTINES
//-------------------------------

double fact(int q);
double binom(int q, int i);
double Inm(int n,int m);

double Inm(int n,int m)
{
  int k,k3,jmin,jmax,j;
  double sum;

  k  = n-m;
  k3 = int(k/3);

  if ( (k%3) != 0 ){
    return 0.;
  }   
  else{
     sum  = 0;
     jmin = 0;
     if ( k3>0 ) jmin = k3;
     jmax = n/3;
     
     for (j=jmin; j<=jmax; j++){

        sum += 2*fact(n)*fact(m)*binom(3*(n-j-k3+1),n-3*j)/
               (fact(n-j-k3+1)*fact(n-j-k3+2)*fact(j-k3)*fact(j) );
     }
  }
 
  return sum;
}
//____________________________________________________________________

double fact(int q)
{
  int    j;
  double f;

  f = 1.0;
  for (j = 2; j<=q; j++) f *= j;

  return f;

}
//____________________________________________________________________

double binom(int q,int i)
{
  int    j;
  double num,den;

  num = 1.0;
  den = 1.0;

  for (j = 1; j<=i; j++){
    num *= (q+1-j);
    den *= j;
  }

  return (num/den);

}

#endif 
