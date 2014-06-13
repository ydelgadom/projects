#ifndef _LOCAL_UPDATE_H
#define _LOCAL_UPDATE_H

#define ERROR_PLAQ(x) \
    if (abs((x))>=LENGTH) { cout << "FATAL ERROR plaq: plaq>max_LENGTH" << endl; return -1; }

#define ERROR_LINK(x,y) \
          if ( ((x)>=LENGTH) || ((y)>=LENGTH) ){ \
            cout << "FATAL ERROR sweep lprime: flux0 > max_LENGTH " << (x) << " " << (y) << endl; \
            return -1; \
          }

//_________________________________________________________________________
int sweep_lprime( int lnsite, int neib[][8]  )
{
  int is0, nu;
  for ( is0=0; is0<lnsite; ++is0 )
  {
    for ( nu=0; nu<4; ++nu )
    {
      int ifield;
      for (ifield=0; ifield<2; ++ifield )
      {
        double ran[2];
        ranlxd(ran,2);

        int ll = vllink[ is0 ][ nu + 4*ifield ];
        int delta = 1-2*(ran[0]<0.5);
        int lnew = ll + delta;
      
        if( lnew>=0 )
        {
          int is1 = neib[is0][nu];
          int ff0 = vflux[is0][ifield];
          int ff1 = vflux[is1][ifield];
          int ffnew0 = ff0 + 2*delta;
          int ffnew1 = ff1 + 2*delta;
          int absk = abs( vlink[is0][nu+4*ifield] );

					ERROR_LINK(ffnew0,ffnew1);

          double weig =   Pn[ffnew0] - Pn[ff0]
                        + Pn[ffnew1] - Pn[ff1]
                        - fac[lnew]  + fac[ll]
                        - fac[absk+lnew]  + fac[absk+ll];
          weig = exp( weig );//*/

          if( ran[1] <= weig )
          {
            vllink[ is0 ][ nu + 4*ifield ] = lnew;
            vflux[ is0 ][ ifield ] = ffnew0;
            vflux[ is1 ][ ifield ] = ffnew1;

            --nflux[ifield][ ff0 ];
            --nflux[ifield][ ff1 ];

            ++nflux[ifield][ ffnew0 ];
            ++nflux[ifield][ ffnew1 ];
          }
        }

      } // field
    } // nu
  } // is

  return 1;
}

#ifdef CUBES
//_________________________________________________________________________
inline void update_a(int *a, int* im, int *ip)
{
	a[abs(im[0])]--;
	a[abs(im[1])]--;
	a[abs(im[2])]--;
	a[abs(im[3])]--;
	a[abs(im[4])]--;
	a[abs(im[5])]--;

	a[abs(ip[0])]++;
	a[abs(ip[1])]++;
	a[abs(ip[2])]++;
	a[abs(ip[3])]++;
	a[abs(ip[4])]++;
	a[abs(ip[5])]++;
}

//_________________________________________________________________________
inline double weight(int *ppnew, int *pp )
{
	double weig = In[abs(ppnew[0])] - In[abs(pp[0])] + 
			          In[abs(ppnew[1])] - In[abs(pp[1])] + 
			          In[abs(ppnew[2])] - In[abs(pp[2])] + 
			          In[abs(ppnew[3])] - In[abs(pp[3])] + 
			          In[abs(ppnew[4])] - In[abs(pp[4])] + 
			          In[abs(ppnew[5])] - In[abs(pp[5])];
   weig = exp( weig );

	return weig;
}

//_________________________________________________________________________
void sweep_cube012( )
{
  int is;
  for(is=0; is<nsite; is++)
  {
    double ran[2];
    ranlxd(ran,2);
    int delta = 1-2*(ran[0]<0.5);

    int xx[4];
    xx[0] = is;
    xx[1] = neib[is][0];
    xx[2] = neib[is][1];
    xx[3] = neib[is][2];

    int pp[6];
    pp[0] = vplaq[xx[0]][0+2-1];
    pp[1] = vplaq[xx[1]][1+2-1];
    pp[2] = vplaq[xx[2]][0+2-1];
    pp[3] = vplaq[xx[0]][1+2-1];
    pp[4] = vplaq[xx[3]][0+1-1];
    pp[5] = vplaq[xx[0]][0+1-1];

    int ppnew[6];
    ppnew[0] = pp[0] - delta; ERROR_PLAQ(ppnew[0]);
    ppnew[1] = pp[1] - delta; ERROR_PLAQ(ppnew[1]);
    ppnew[2] = pp[2] + delta; ERROR_PLAQ(ppnew[2]);
    ppnew[3] = pp[3] + delta; ERROR_PLAQ(ppnew[3]);
    ppnew[4] = pp[4] - delta; ERROR_PLAQ(ppnew[4]);
    ppnew[5] = pp[5] + delta; ERROR_PLAQ(ppnew[5]);

    if( ran[1]<=weight( ppnew, pp ) )
    {
      vplaq[xx[0]][0+2-1] = ppnew[0];
      vplaq[xx[1]][1+2-1] = ppnew[1];
      vplaq[xx[2]][0+2-1] = ppnew[2];
      vplaq[xx[0]][1+2-1] = ppnew[3];
      vplaq[xx[3]][0+1-1] = ppnew[4];
      vplaq[xx[0]][0+1-1] = ppnew[5];

			update_a(nplaq, pp, ppnew);
    }
  }

}

//_________________________________________________________________________

void sweep_cubexy3( int x, int y )
{
  int is;
  for(is=0; is<nsite; is++)
  {
    double ran[2];
    ranlxd(ran,2);
    int delta = 1-2*(ran[0]<0.5);

    int xx[4], pp[6], ppnew[6];
    xx[0] = is;
    xx[1] = neib[is][x];
    xx[2] = neib[is][y];
    xx[3] = neib[is][3];

    pp[0] = vplaq[xx[0]][x+3];
    pp[1] = vplaq[xx[1]][y+3];
    pp[2] = vplaq[xx[2]][x+3];
    pp[3] = vplaq[xx[0]][y+3];
    pp[4] = vplaq[xx[3]][x+y-1];
    pp[5] = vplaq[xx[0]][x+y-1];

    ppnew[0] = pp[0] - delta; ERROR_PLAQ(ppnew[0]);
    ppnew[1] = pp[1] - delta; ERROR_PLAQ(ppnew[1]);
    ppnew[2] = pp[2] + delta; ERROR_PLAQ(ppnew[2]);
    ppnew[3] = pp[3] + delta; ERROR_PLAQ(ppnew[3]);
    ppnew[4] = pp[4] - delta; ERROR_PLAQ(ppnew[4]);
    ppnew[5] = pp[5] + delta; ERROR_PLAQ(ppnew[5]);
 
    if( ran[1]<=weight( ppnew, pp ) )
    {
      vplaq[xx[0]][x+3] = ppnew[0];
      vplaq[xx[1]][y+3] = ppnew[1];
      vplaq[xx[2]][x+3] = ppnew[2];
      vplaq[xx[0]][y+3] = ppnew[3];
      vplaq[xx[3]][x+y-1] = ppnew[4];
      vplaq[xx[0]][x+y-1] = ppnew[5];

			update_a(nplaq, pp, ppnew);
    }
  }
}

//________________________________________________________________________
void sweep_cube( )
{
  sweep_cube012( );
  sweep_cubexy3( 0, 1 );
	sweep_cubexy3( 0, 2 );
	sweep_cubexy3( 1, 2 );
}

#endif


#ifdef WL
//_________________________________________________________________________
int sweep_winding_loop( int vleng, int neib[][8] )
{
  int leng2 = vleng*vleng;
  int leng3 = vleng*leng2;

  for ( int i3=0; i3<leng_t; i3++){
    for ( int i2=0; i2<vleng; i2++ ){
      for ( int i1=0; i1<vleng; i1++ ){
        double ran[2];
        ranlxd(ran,2);
        int ix = i3*leng3 + i2*leng2 + i1*vleng;
        int delta = 1-2*(ran[0]<0.5);

        int link[vleng][2];
        link[vleng-1][0] = vlink[ix + vleng-1][0];
        link[vleng-1][1] = vlink[ix + vleng-1][4];
        int i0; double weig = 0.;
        for ( i0=0; i0<vleng; i0++ )
        {
          int xx = ix + i0;
          int iold = (i0-1)*(i0!=0) + (i0==0)*(vleng-1);

          link[i0][0] = vlink[xx][0];
          link[i0][1] = vlink[xx][4];

          int ffnew[2];
          ffnew[0] = vflux[xx][0] + abs(link[i0][0] + delta) - abs(link[i0][0]) + abs(link[iold][0] + delta) - abs(link[iold][0]);
          ffnew[1] = vflux[xx][1] + abs(link[i0][1] + delta) - abs(link[i0][1]) + abs(link[iold][1] + delta) - abs(link[iold][1]);

					ERROR_LINK(ffnew[0],ffnew[1]);

          weig += ( Pn[ffnew[0]] - Pn[vflux[xx][0]] + Pn[ffnew[1]] - Pn[vflux[xx][1]]
                    + fac[ abs(link[i0][0]) + vllink[xx][0] ] - fac[ abs(link[i0][0] + delta) + vllink[xx][0] ]
                    + fac[ abs(link[i0][1]) + vllink[xx][4] ] - fac[ abs(link[i0][1] + delta) + vllink[xx][4] ] );
        }
        weig = exp( weig );

        if ( ran[1] < weig ){
          for ( i0=0; i0<vleng; i0++ )
          {
            int xx = ix + i0;
            int iold = (i0-1)*(i0!=0) + (i0==0)*(vleng-1);

            vlink[ xx ][ 0 ] += delta;
            vlink[ xx ][ 4 ] += delta;

            nflux[0][ vflux[xx][0] ]--;
            nflux[1][ vflux[xx][1] ]--;
            vflux[ xx ][ 0 ] += (abs(link[i0][0] + delta) - abs(link[i0][0]) + abs(link[iold][0]+delta) - abs(link[iold][0]));
            vflux[ xx ][ 1 ] += (abs(link[i0][1] + delta) - abs(link[i0][1]) + abs(link[iold][1]+delta) - abs(link[iold][1]));
            nflux[0][ vflux[xx][0] ]++;
            nflux[1][ vflux[xx][1] ]++;
          }
        }
      }
    }
  }

  ////////////////////////
  for ( int i3=0; i3<leng_t; i3++) {
    for ( int i2=0; i2<vleng; i2++ ){
      for ( int i0=0; i0<vleng; i0++ ){
        double ran[2];
        ranlxd(ran,2);
        int ix = i3*leng3 + i2*leng2 + i0;
        int delta = (1-2*(ran[0]<0.5));

        int link[vleng][2];
        link[vleng-1][0] = vlink[ix + leng2-vleng][1];
        link[vleng-1][1] = vlink[ix + leng2-vleng][5];
        int i1; double weig = 0.;
        for ( i1=0; i1<vleng; i1++ )
        {
          int xx = ix + i1*vleng;
          int iold = (i1-1)*(i1!=0) + (i1==0)*(vleng-1);

          link[i1][0] = vlink[xx][1];
          link[i1][1] = vlink[xx][5];

          int ffnew[2];
          ffnew[0] = vflux[xx][0] + abs(link[i1][0] + delta) - abs(link[i1][0]) + abs(link[iold][0] + delta) - abs(link[iold][0]);
          ffnew[1] = vflux[xx][1] + abs(link[i1][1] + delta) - abs(link[i1][1]) + abs(link[iold][1] + delta) - abs(link[iold][1]);

          ERROR_LINK(ffnew[0],ffnew[1]);

          weig += ( Pn[ffnew[0]] - Pn[vflux[xx][0]] + Pn[ffnew[1]] - Pn[vflux[xx][1]]
                    + fac[ abs(link[i1][0]) + vllink[xx][1] ] - fac[ abs(link[i1][0] + delta) + vllink[xx][1] ]
                    + fac[ abs(link[i1][1]) + vllink[xx][5] ] - fac[ abs(link[i1][1] + delta) + vllink[xx][5] ] );
        }
        weig = exp( weig );

        if ( ran[1] < weig ){
          for ( i1=0; i1<vleng; i1++ )
          {
            int xx = ix + i1*vleng;
            int iold = (i1-1)*(i1!=0) + (i1==0)*(vleng-1);

            vlink[ xx ][ 1 ] += delta;
            vlink[ xx ][ 5 ] += delta;

            nflux[0][ vflux[xx][0] ]--;
            nflux[1][ vflux[xx][1] ]--;
            vflux[ xx ][ 0 ] += (abs(link[i1][0] + delta) - abs(link[i1][0]) + abs(link[iold][0]+delta) - abs(link[iold][0]));
            vflux[ xx ][ 1 ] += (abs(link[i1][1] + delta) - abs(link[i1][1]) + abs(link[iold][1]+delta) - abs(link[iold][1]));
            nflux[0][ vflux[xx][0] ]++;
            nflux[1][ vflux[xx][1] ]++;
          }
        }
      }
    }
  }

  /////
  for ( int i3=0; i3<leng_t; i3++){
    for ( int i1=0; i1<vleng; i1++ ){
      for ( int i0=0; i0<vleng; i0++ ){
        double ran[2];
        ranlxd(ran,2);
        int ix = i3*leng3 + i1*vleng + i0;
        int delta = 1-2*(ran[0]<0.5);

        int link[vleng][2];
        link[vleng-1][0] = vlink[ix + leng3-leng2][2];
        link[vleng-1][1] = vlink[ix + leng3-leng2][6];
        int i2; double weig = 0.;
        for ( i2=0; i2<vleng; i2++ )
        {
          int xx = ix + i2*leng2;
          int iold = (i2-1)*(i2!=0) + (i2==0)*(vleng-1);

          link[i2][0] = vlink[xx][2];
          link[i2][1] = vlink[xx][6];

          int ffnew[2];
          ffnew[0] = vflux[xx][0] + abs(link[i2][0] + delta) - abs(link[i2][0]) + abs(link[iold][0] + delta) - abs(link[iold][0]);
          ffnew[1] = vflux[xx][1] + abs(link[i2][1] + delta) - abs(link[i2][1]) + abs(link[iold][1] + delta) - abs(link[iold][1]);

					ERROR_LINK(ffnew[0],ffnew[1]);

          weig += ( Pn[ffnew[0]] - Pn[vflux[xx][0]] + Pn[ffnew[1]] - Pn[vflux[xx][1]]
                    + fac[ abs(link[i2][0]) + vllink[xx][2] ] - fac[ abs(link[i2][0] + delta) + vllink[xx][2] ]
                    + fac[ abs(link[i2][1]) + vllink[xx][6] ] - fac[ abs(link[i2][1] + delta) + vllink[xx][6] ] );
        }
        weig = exp( weig );

        if ( ran[1] < weig ){
          for ( i2=0; i2<vleng; i2++ )
          {
            int xx = ix + i2*leng2;
            int iold = (i2-1)*(i2!=0) + (i2==0)*(vleng-1);

            vlink[ xx ][ 2 ] += delta;
            vlink[ xx ][ 6 ] += delta;

            nflux[0][ vflux[xx][0] ]--;
            nflux[1][ vflux[xx][1] ]--;
            vflux[ xx ][ 0 ] += (abs(link[i2][0] + delta) - abs(link[i2][0]) + abs(link[iold][0]+delta) - abs(link[iold][0]));
            vflux[ xx ][ 1 ] += (abs(link[i2][1] + delta) - abs(link[i2][1]) + abs(link[iold][1]+delta) - abs(link[iold][1]));
            nflux[0][ vflux[xx][0] ]++;
            nflux[1][ vflux[xx][1] ]++;
          }
        }
      }
    }
  }

  ///////
  for ( int i2=0; i2<vleng; i2++ ){
    for ( int i1=0; i1<vleng; i1++ ){
      for ( int i0=0; i0<vleng; i0++){
        double ran[2];
        ranlxd(ran,2);
        int ix = i2*leng2 + i1*vleng + i0;
        int delta = (1-2*(ran[0]<0.5));

        int link[leng_t][2];
        link[leng_t-1][0] = vlink[ix + leng3*(leng_t-1)][3];
        link[leng_t-1][1] = vlink[ix + leng3*(leng_t-1)][7];
        double weig = 0.;
        for ( int i3=0; i3<leng_t; i3++ )
        {
          int xx = ix + leng3*i3;
          int iold = (i3-1)*(i3!=0) + (i3==0)*(leng_t-1);

          link[i3][0] = vlink[xx][3];
          link[i3][1] = vlink[xx][7];

          int ffnew[2];
          ffnew[0] = vflux[xx][0] + abs(link[i3][0] + delta) - abs(link[i3][0]) + abs(link[iold][0] + delta) - abs(link[iold][0]);
          ffnew[1] = vflux[xx][1] + abs(link[i3][1] + delta) - abs(link[i3][1]) + abs(link[iold][1] + delta) - abs(link[iold][1]);

					ERROR_LINK(ffnew[0],ffnew[1]);

          weig += ( Pn[ffnew[0]] - Pn[vflux[xx][0]] + Pn[ffnew[1]] - Pn[vflux[xx][1]]
                    + fac[ abs(link[i3][0]) + vllink[xx][3] ] - fac[ abs(link[i3][0] + delta) + vllink[xx][3] ]
                    + fac[ abs(link[i3][1]) + vllink[xx][7] ] - fac[ abs(link[i3][1] + delta) + vllink[xx][7] ]
                    + (2*mu)*delta );
        }
        weig = exp( weig );

        if ( ran[1] < weig ){
          for ( int i3=0; i3<leng_t; i3++ ){
            int xx = ix + i3*leng3;
            int iold = (i3-1)*(i3!=0) + (i3==0)*(leng_t-1);

            vlink[ xx ][ 3 ] += delta;
            vlink[ xx ][ 7 ] += delta;
            nlink[0] += delta;
            nlink[1] += delta;

            nflux[0][ vflux[xx][0] ]--;
            nflux[1][ vflux[xx][1] ]--;
            vflux[ xx ][ 0 ] += (abs(link[i3][0] + delta) - abs(link[i3][0]) + abs(link[iold][0]+delta) - abs(link[iold][0]));
            vflux[ xx ][ 1 ] += (abs(link[i3][1] + delta) - abs(link[i3][1]) + abs(link[iold][1]+delta) - abs(link[iold][1]));
            nflux[0][ vflux[xx][0] ]++;
            nflux[1][ vflux[xx][1] ]++;
          }
        }
      }
    }
  }

  return 1;

}

#endif

#endif
