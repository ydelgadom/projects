#ifndef _SWEEPS_H
#define _SWEEPS_H

#define ERROR_SWEEP(x) if ((x)>=LENGTH_FLUX) \
					{ cout << "FATAL ERROR sweep: flux>max_length" << endl; exit(-1); }

void sweep_l( )
{
	/*
		Local update of unconstrained variable
	*/
  for ( int is0=0; is0<nsite; is0++ )
    for ( int nu=0; nu<4; nu++ )
    {
      double ran[2];
      ranlxd(ran,2);

      int ll = vllink[is0][nu];
      int delta = 1-2*(ran[0]<0.5);
      int lnew = ll + delta;
      
      if( lnew>=0 )
      {
        int is1 = neib[is0][nu];
        int ff0 = vflux[is0];
        int ff1 = vflux[is1];
        int ffnew0 = ff0 + 2*delta;
        int ffnew1 = ff1 + 2*delta;
        int absk = abs(vlink[is0][nu]);

				// Check if link variable > maximum
				ERROR_SWEEP(ffnew0);
				ERROR_SWEEP(ffnew1);

        double weig = Pn[ffnew0] - Pn[ff0]
                   + Pn[ffnew1] - Pn[ff1]
                   - fac[lnew]  + fac[ll]
                   - fac[absk+lnew]  + fac[absk+ll];
        weig = exp( weig );//*/

        if( ran[1] <= weig )
        {
          vllink[is0][nu] = lnew;
          vflux[is0] = ffnew0;
          vflux[is1] = ffnew1;

          nflux[ff0]--;
          nflux[ff1]--;

          nflux[ffnew0]++;
          nflux[ffnew1]++;
        }
      }
    }
}

#endif
