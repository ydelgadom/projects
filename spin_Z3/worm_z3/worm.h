#ifndef _WORM_H
#define _WORM_H

/**************************************************************************

  This file contains:

  Close worm:
  void nworms(int nw)
  void add_link()

  Open worm:
  void nworms_open(int nw)
  void add_link_open()

  The pseudo-code of both algorithms is in the CPC papers in ../../papers

***************************************************************************/

int pp,p,x,x0,wormsgn;

#ifndef OPEN
//_________________________________________________________________________
void add_link()
{
  int    xnew,nu,nuadd,xadd;
  int    d,dnew,delta;
  double ran[2];

  if (p){

    do{
      ranlxd(ran,2);
      xnew = int(ran[0]*nsite);
      d = mon[xnew];
      delta = 2 - wormsgn;
      dnew = triadd[d][delta];

      if (ran[1] <= monoweight[dnew]/monoweight[d]){
        x = xnew;
        p = false;
        mon[x] = dnew;
        nmon[d]--;
        nmon[dnew]++;
      }
    }while(p);

  }
  else{
    ranlxd(ran,1);
    nu = int(ran[0]*9);

    if ( nu >= 6 ){
      d = mon[x];
      dnew = triadd[d][wormsgn];
      ranlxd(ran,1);
      if ( ran[0] <= monoweight[dnew]/monoweight[d]){
        p = true;
        mon[x] = dnew;
        nmon[d]--;
        nmon[dnew]++;
      }
    }
    else{

      xnew = neib[x][nu];

      if ( nu < 3 ){
        xadd = x;
        nuadd = nu;
        delta = wormsgn;
      }
      else{
        xadd = xnew;
        nuadd = nu - 3;
        delta = 2-wormsgn;
      }

      d = dim[xadd][nuadd];

      if ( d == 1 ){
        ranlxd(ran,1);
        if ( ran[0] <= bb){
          x = xnew;
          dnew = triadd[d][delta];
          dim[xadd][nuadd] = dnew;
          ndim[d]--;
          ndim[dnew]++;
        }
      }
      else{
        x = xnew;
        dnew = triadd[d][delta];
        dim[xadd][nuadd] = dnew;
        ndim[d]--;
        ndim[dnew]++;
      }
    }
  }

}

//_________________________________________________________________________
void nworms(int nw)
{
  int   iw;
  double ran2[2];

  for(iw=1; iw<=nw; ++iw){
    ranlxd(ran2,2);
    
    wormsgn = 2; // increase 1
    if ( ran2[0] < 0.5 ) wormsgn = 0;  // decrease 1
   
    p  = false;
    x  = int(ran2[1]*nsite);
    x0 = x;

    do{
      add_link();
    }while ( (x != x0) || p );
  }

}
#endif

#ifdef OPEN
//_________________________________________________________________________
void add_link_open()
{
  int    xnew,nu,nuadd,xadd;
  int    d,dnew,delta;
  double ran[2];

  ranlxd(ran,1);
  if (pp==0) 
    nu = 6;
  else
    nu = int(ran[0]*7);

  if ( nu > 5 )
  {
    if (x==x0 && pp==0)
      delta = 2 - wormsgn;
    else
      delta = wormsgn;

    d    = mon[x];
    dnew = triadd[d][delta];
    ranlxd(ran,1);
    if ( ran[0] <= monoweight[dnew]/monoweight[d])
    {
      pp++;
      mon[x] = dnew;
      nmon[d]--;
      nmon[dnew]++;
    }
    if ( pp==0 && x==x0 ) pp = 2;
  }
  else
  {
    xnew = neib[x][nu];

    if ( nu < 3 ){
      xadd = x;
      nuadd = nu;
      delta = wormsgn;
    }
    else{
      xadd = xnew;
      nuadd = nu - 3;
      delta = 2-wormsgn;
    }

    d = dim[xadd][nuadd];
    if ( d == 1 )
    {
      ranlxd(ran,1);
      if ( ran[0] <= bb){
        x = xnew;
        dnew = triadd[d][delta];
        dim[xadd][nuadd] = dnew;
        ndim[d]--;
        ndim[dnew]++;
      }
    }
    else{
      x = xnew;
      dnew = triadd[d][delta];
      dim[xadd][nuadd] = dnew;
      ndim[d]--;
      ndim[dnew]++;
    }
  }

}

//_________________________________________________________________________
void nworms(int nw)
{
  int   iw;
  double ran2[2];

  for(iw=1; iw<=nw; ++iw){
  
    ranlxd(ran2,2);
    
    wormsgn = 2; // increase 1
    if ( ran2[0] < 0.5 ) wormsgn = 0;  // decrease 1
   
    pp = 0;
    x  = int(ran2[1]*nsite);
    x0 = x;

    do{
      add_link_open();
    }while ( pp!=2 );
  }

}
#endif

#endif
