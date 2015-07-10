#ifndef _WORM_H
#define _WORM_H

/************************************************************************* 
*   This file contains:
*
*   Total update:
*   void nworms(int nw)
*
*   Worm update for constrained variables:
*   void move0()
*   void move1()
*   void move()
*   void get_current_segments_and_links( int rho, int sseg, int *xf )
*
*   Update of unconstrained variables:
*   void sweep_l( )
*
*------------------------------------------------------------------------------
*
* Directions:
* nu,rho = 0,1,2 spatial directions
* nu,rho = 3     temporal direction
*
* Link ordering in the plaquettes:
*         __2__ 
*        |     |
*     ^  3     1
*  nu1|  |__0__| 
*     --->
*     nu0     nu0<nu1 always
*
* vlink[i][j]: i=site and j=0,1,2,3
* vplaq[i][j]: i=site and j=nu+rho-1 if spatial plaquette
*              i=site and j=nu+rho   if temporal plaquette
*
* xx  = site of open_link_init
* nu  = direction of open_link_init :   nu>0  always
* rho = other direction of plaquette:   rho>0 moving towards positive direction of axis
*                                       rho<0 moving towards negative direction of axis
*
* open_link_init: (x) ---> (x+nu)
* open_link_end:  choose randomly among the other 3 links
*
* Regarding variables:
* ivar = index of var(=segment or link)
* vvar = value of var(=segment or link)
*
**************************************************************************/

#define NWORMS nsite

#define ERROR_WORM(x) \
      if ((x)>=LENGTH_FLUX) \
          { cout << "FATAL ERROR WORM: flux>max_length" << endl; exit(-1); }

#define ERROR_PLAQ(x) \
    if (abs((x))>=LENGTH) { cout << "FATAL ERROR plaq: plaq>max_length" << endl; exit(-1); }

// STRUCT VARS
struct tlink
{
  /* Link variable */
  int lidx;// position in the segment [0,3]
  int lval; // occupation number
  int lx;   // lattice site [0,nsites-1]
  int lnu;  // link direction [0,3]
};

struct segment
{
  /* Segment variable for worm updates */
  tlink slink[4]; // 4 link variables: link[0] and link[1] = links to be modified 
                  //                   link[2] = current worm's head 
                  //                   link[3] = new worm's head
  int   svplaq;   // occupation number of plaquette variable
  int   sxplaq;   // site of the lower-left corner of the plaquette
  int   sdirmove; // if the segment will be inserted in positive or negative direction
  int   ssignseg; // if a positve or negative segment will be inserted
};

// Global variables 
int    insert_next_segment;
int    oldsignseg,oldrho,real_world,oldlidx2;
int    delta_seg[5] = {-1,-1,1,1,1}; //  link_0,link_1,link_2,link_3,pp = -1,-1,+1,+1,+1
tlink   link0;
segment seg;

// Subroutines
//_________________________________________________________________________
inline void a_eq_b(int *a, int*b )
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
  a[3] = b[3];
}

//________________________________________________________________________
inline void a_eq_mb(int *a, int*b )
{
  a[0] = -b[0];
  a[1] = -b[1];
  a[2] = -b[2];
  a[3] = -b[3];
}

//________________________________________________________________________
void get_current_segments_and_links( int rho, int sseg, int *xf )
{
  // get 4 links and plaquette and set
  // open_link_init and open_link_end

  if (seg.slink[2].lnu<rho)
  {
    if (seg.sdirmove > 0)
    {
      seg.slink[2].lidx = 0;
      seg.slink[3].lidx++;
      switch (seg.slink[3].lidx)
      {
        case 1:
        seg.ssignseg = 0;

        seg.slink[0].lx = neib[seg.slink[2].lx][rho];
        seg.slink[0].lnu = seg.slink[2].lnu;
        seg.slink[0].lidx = 2;

        seg.slink[1].lx = seg.slink[2].lx;
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 3;

        seg.slink[3].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[3].lnu = rho;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = -1;
        break;
        
        case 2:
        seg.ssignseg = 0;

        seg.slink[0].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 1;

        seg.slink[1].lx = seg.slink[2].lx;
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 3;

        seg.slink[3].lx = neib[seg.slink[2].lx][rho];
        seg.slink[3].lnu = seg.slink[2].lnu;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        break;

        case 3:
        seg.ssignseg = 1;

        seg.slink[0].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 1;

        seg.slink[1].lx = neib[seg.slink[2].lx][rho];
        seg.slink[1].lnu = seg.slink[2].lnu;
        seg.slink[1].lidx = 2;

        seg.slink[3].lx = seg.slink[2].lx;
        seg.slink[3].lnu = rho;

        xf[0] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[1] = seg.slink[0].lx;
        xf[2] = seg.slink[1].lx;
        xf[3] = -1;
        break;        
      }
      seg.sxplaq = seg.slink[2].lx;
    }
    else
    {
      seg.slink[2].lidx = 2;
      seg.slink[3].lidx = seg.slink[2].lidx - seg.slink[3].lidx + (seg.slink[3].lidx==0);
      switch (seg.slink[3].lidx)
      {
        case 0:
        seg.ssignseg = 1;

        seg.slink[3].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[3].lnu = seg.slink[2].lnu;
        seg.sxplaq = seg.slink[3].lx;

        seg.slink[0].lx = neib[seg.slink[3].lx][seg.slink[2].lnu];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 1;

        seg.slink[1].lx = seg.slink[3].lx;
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 3;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        break;
        
        case 1:
        seg.ssignseg = 0;

        seg.slink[0].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[0].lnu = seg.slink[2].lnu;
        seg.slink[0].lidx = 0;
        seg.sxplaq = seg.slink[0].lx;

        seg.slink[1].lx = seg.slink[0].lx;
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 3;

        seg.slink[3].lx = neib[seg.slink[0].lx][seg.slink[2].lnu];
        seg.slink[3].lnu = rho;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[3] = -1;
        break;

        case 3:
        seg.ssignseg = 1;

        seg.slink[3].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[3].lnu = rho;

        seg.slink[0].lx = seg.slink[3].lx;
        seg.slink[0].lnu = seg.slink[2].lnu;
        seg.slink[0].lidx = 0;
        seg.sxplaq = seg.slink[0].lx;

        seg.slink[1].lx = neib[seg.slink[3].lx][seg.slink[2].lnu];
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 1;

        xf[0] = seg.slink[1].lx;
        xf[1] = seg.slink[0].lx;
        xf[2] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[3] = -1;
        break;        
      }
    }
  }
  else
  {
    if (seg.sdirmove > 0)
    {
      seg.slink[2].lidx = 3;
      switch (seg.slink[3].lidx)
      {
        case 0:
        seg.ssignseg = 1;

        seg.slink[3].lx = seg.slink[2].lx;
        seg.slink[3].lnu = rho;

        seg.slink[0].lx = neib[seg.slink[2].lx][rho];
        seg.slink[0].lnu = seg.slink[2].lnu;
        seg.slink[0].lidx = 1;

        seg.slink[1].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 2;

        xf[0] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[1] = seg.slink[0].lx;
        xf[2] = seg.slink[1].lx;
        xf[3] = -1;
        break;
        
        case 1:
        seg.ssignseg = 0;

        seg.slink[0].lx = seg.slink[2].lx;
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 0;

        seg.slink[1].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 2;

        seg.slink[3].lx = neib[seg.slink[2].lx][rho];
        seg.slink[3].lnu = seg.slink[2].lnu;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        break;

        case 2:
        seg.ssignseg = 0;

        seg.slink[0].lx = seg.slink[2].lx;
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 0;

        seg.slink[1].lx = neib[seg.slink[2].lx][rho];
        seg.slink[1].lnu = seg.slink[2].lnu;
        seg.slink[1].lidx = 1;

        seg.slink[3].lx = neib[seg.slink[2].lx][seg.slink[2].lnu];
        seg.slink[3].lnu = rho;

        xf[0] = seg.slink[1].lx;
        xf[1] = seg.slink[0].lx;
        xf[2] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[3] = -1;
        break;        
      }
      seg.sxplaq = seg.slink[2].lx;     
    }
    else
    {
      seg.slink[2].lidx = 1;
      seg.slink[3].lidx += (seg.slink[3].lidx!=0);
      switch (seg.slink[3].lidx)
      {
        case 0:
        seg.ssignseg = 1;

        seg.slink[3].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[3].lnu = rho;
        seg.sxplaq = seg.slink[3].lx;

        seg.slink[0].lx = neib[seg.slink[3].lx][seg.slink[2].lnu];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 2;

        seg.slink[1].lx = seg.slink[3].lx;
        seg.slink[1].lnu = seg.slink[2].lnu;
        seg.slink[1].lidx = 3;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = -1;
        break;
        
        case 2:
        seg.ssignseg = 0;

        seg.slink[0].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 0;
        seg.sxplaq = seg.slink[0].lx;

        seg.slink[1].lx = seg.slink[0].lx;
        seg.slink[1].lnu = seg.slink[2].lnu;
        seg.slink[1].lidx = 3;

        seg.slink[3].lx = neib[seg.slink[0].lx][seg.slink[2].lnu];
        seg.slink[3].lnu = rho;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        xf[3] = -1;
        break;

        case 3:
        seg.ssignseg = 1;

        seg.slink[0].lx = neib[seg.slink[2].lx][rho+4];
        seg.slink[0].lnu = rho;
        seg.slink[0].lidx = 0;
        seg.sxplaq = seg.slink[0].lx;

        seg.slink[1].lx = neib[seg.slink[0].lx][seg.slink[2].lnu];
        seg.slink[1].lnu = rho;
        seg.slink[1].lidx = 2;

        seg.slink[3].lx = seg.slink[0].lx;
        seg.slink[3].lnu = seg.slink[2].lnu;

        xf[0] = seg.slink[0].lx;
        xf[1] = neib[seg.slink[0].lx][seg.slink[0].lnu];
        xf[2] = seg.slink[1].lx;
        xf[3] = neib[seg.slink[1].lx][seg.slink[1].lnu];
        break;        
      }      
    }
  }
  seg.slink[0].lval = vlink[seg.slink[0].lx][seg.slink[0].lnu];
  seg.slink[1].lval = vlink[seg.slink[1].lx][seg.slink[1].lnu];
  seg.slink[3].lval = vlink[seg.slink[3].lx][seg.slink[3].lnu];
  seg.svplaq = vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg];
}

//_________________________________________________________________________
void move0()
{
  //
  // Insert first link
  //
  double ran[4];
  ranlxd( ran, 4 );
  seg.slink[2].lx = int(ran[0]*nsite);
  seg.slink[2].lnu = int(ran[1]*4);
  seg.slink[2].lval = vlink[seg.slink[2].lx][seg.slink[2].lnu];

  // calculate probability
  link0.lval = 1-2*(ran[2]<0.5);
  int newvalseg = seg.slink[2].lval + link0.lval;

  int ll = vllink[seg.slink[2].lx][seg.slink[2].lnu];
  int ff[2], ffnew[2];
  ff[0] = vflux[seg.slink[2].lx];
  ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ];

  ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
  ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);

  // Check if fluxes > maximum
  ERROR_WORM(ffnew[0]);
  ERROR_WORM(ffnew[1]);

  double weig = Pn[ffnew[0]] - Pn[ff[0]]
         + Pn[ffnew[1]] - Pn[ff[1]]
         + fac[ abs(seg.slink[2].lval) + ll ] - fac[ abs(newvalseg) + ll ];
  weig = exp( weig );

  if ( ran[3] <= weig )
  {
    insert_next_segment = 1;
    
    vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg;
    vflux[seg.slink[2].lx] = ffnew[0];
    vflux[neib[seg.slink[2].lx][seg.slink[2].lnu]] = ffnew[1];

    nflux[ff[0]]--;
    nflux[ff[1]]--;

    nflux[ffnew[0]]++;
    nflux[ffnew[1]]++;

    nlink[abs(seg.slink[2].lval)]--;
    nlink[abs(newvalseg)]++;
  }

}

//_________________________________________________________________________
void move1()
{
  // 
  // The worm tries to insert the first segment
  //
  double ran[2];
  ranlxd(ran,2);
  int rho = int(ran[0]*4); 

  if ( rho == 3 ) // try to close the worm
  {
    seg.slink[2].lval = vlink[seg.slink[2].lx][seg.slink[2].lnu];
    int newvalseg = seg.slink[2].lval - link0.lval;

    int ll = vllink[seg.slink[2].lx][seg.slink[2].lnu];
    int ff[2], ffnew[2];
    ff[0] = vflux[seg.slink[2].lx];
    ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ];

    ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
    ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);

    // Check if fluxes > maximum
    ERROR_WORM(ffnew[0]);
    ERROR_WORM(ffnew[1]);

    double weig = Pn[ffnew[0]] - Pn[ff[0]]
           + Pn[ffnew[1]] - Pn[ff[1]]
           + fac[ abs(seg.slink[2].lval) + ll ] - fac[ abs(newvalseg) + ll ];
    weig = exp( weig );

    if ( ran[1] <= weig )
    {
      vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg;
      vflux[seg.slink[2].lx] = ffnew[0];
      vflux[neib[seg.slink[2].lx][seg.slink[2].lnu] ] = ffnew[1];

      nlink[ abs(seg.slink[2].lval) ]--;
      nlink[ abs(newvalseg) ]++;

      nflux[ff[0]]--;
      nflux[ff[1]]--;

      nflux[ffnew[0]]++;
      nflux[ffnew[1]]++;
        
      insert_next_segment = 0;
    }

  } else  // Insert segment
  {
    seg.sdirmove = ran[1]<0.5; // 1=positive and 0=negative direction
    switch(seg.slink[2].lnu){
      case 0:
      rho++;
      break;

      case 1:
      rho += (rho!=0);
      break;

      case 2:
      rho = seg.slink[2].lnu - rho + (rho==0);
      break;
    }
    ranlxd(ran,2);
    seg.slink[3].lidx = int(ran[0]*3);
    int sseg = (seg.slink[2].lnu!=3)&(rho!=3);
    int xf[4];
    get_current_segments_and_links( rho, sseg, xf );

    // select sign:
    int delta[5];
    if ( link0.lval != delta_seg[seg.slink[2].lidx] )
      a_eq_mb( delta, delta_seg );
    else
      a_eq_b( delta, delta_seg );

    int newvalseg[3];
    newvalseg[0] = seg.slink[0].lval + delta[seg.slink[0].lidx];
    newvalseg[1] = seg.slink[1].lval + delta[seg.slink[1].lidx];
    newvalseg[2] = seg.svplaq + delta[4];

    int ll[2], ff[4], ffnew[4];
    ll[0] = vllink[seg.slink[0].lx][seg.slink[0].lnu];
    ll[1] = vllink[seg.slink[1].lx][seg.slink[1].lnu];

    if ( xf[3]<0 )
    {
      ff[0] = vflux[xf[0]];
      ff[1] = vflux[xf[1]];
      ff[2] = vflux[xf[2]];
      ff[3] = 0;

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval)
                       + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = 0;
    }
    else
    {
      ff[0] = vflux[xf[0]];
      ff[1] = vflux[xf[1]];
      ff[2] = vflux[xf[2]];
      ff[3] = vflux[xf[3]];

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = ff[3] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
    }
    // Check if fluxes > maximum
    ERROR_PLAQ(newvalseg[2]);
    ERROR_WORM(ffnew[0]); 
    ERROR_WORM(ffnew[1]);
    ERROR_WORM(ffnew[2]);
    ERROR_WORM(ffnew[3]);

    // calculate probability
    double weig = In[abs(newvalseg[2])]  - In[abs(seg.svplaq)]
               + Pn[ffnew[0]] - Pn[ff[0]]
               + Pn[ffnew[1]] - Pn[ff[1]]
               + Pn[ffnew[2]] - Pn[ff[2]]
               + Pn[ffnew[3]] - Pn[ff[3]]
               + fac[abs(seg.slink[0].lval)+ll[0]] - fac[abs(newvalseg[0])+ll[0]]
               + fac[abs(seg.slink[1].lval)+ll[1]] - fac[abs(newvalseg[1])+ll[1]];
    weig = exp( weig );

    if ( ran[1] <= weig )
    {
      insert_next_segment = 2;

      vlink[seg.slink[0].lx][seg.slink[0].lnu] = newvalseg[0];
      vlink[seg.slink[1].lx][seg.slink[1].lnu] = newvalseg[1];
      vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg] = newvalseg[2];

      vflux[xf[0]] = ffnew[0];
      vflux[xf[1]] = ffnew[1];
      vflux[xf[2]] = ffnew[2];

      nlink[abs(seg.slink[0].lval)]--;
      nlink[abs(seg.slink[1].lval)]--;
      nplaq[abs(seg.svplaq)]--;

      nlink[abs(newvalseg[0])]++;
      nlink[abs(newvalseg[1])]++;
      nplaq[abs(newvalseg[2])]++;

      nflux[ff[0]]--;
      nflux[ff[1]]--;
      nflux[ff[2]]--;

      nflux[ffnew[0]]++;
      nflux[ffnew[1]]++;
      nflux[ffnew[2]]++;

      if ( xf[3]>=0 )
      {
        vflux[xf[3]] = ffnew[3];
        nflux[ff[3]]--;
        nflux[ffnew[3]]++;
      }

      oldsignseg = seg.ssignseg;
      if (seg.slink[3].lnu == rho){
        oldrho = seg.slink[2].lnu;
      }else{
        oldrho = rho;
      }
      // open_end --> open_init of next move
      seg.slink[2].lx = seg.slink[3].lx;
      seg.slink[2].lnu = seg.slink[3].lnu;
      seg.slink[2].lval = seg.slink[3].lval;
      oldlidx2 = seg.slink[3].lidx;

      a_eq_b( delta_seg, delta );
    }// accepted update
  }
}

//_________________________________________________________________________
void move( )
{
  double ran[2];
  ranlxd(ran,2);
  int rho = int(ran[0]*4); 

  if ( rho == 3 ) // try to close worm
  {
    int newvalseg = seg.slink[2].lval + delta_seg[oldlidx2];

    int ll = vllink[seg.slink[2].lx][seg.slink[2].lnu];
    int ff[2], ffnew[2];
    ff[0] = vflux[seg.slink[2].lx];
    ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ];

    ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
    ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);

    ERROR_WORM(ffnew[0]); 
    ERROR_WORM(ffnew[1]);

    double weig = Pn[ffnew[0]] - Pn[ff[0]]
           + Pn[ffnew[1]] - Pn[ff[1]]
           + fac[ abs(seg.slink[2].lval) + ll ] - fac[ abs(newvalseg) + ll ];
    weig = exp( weig );

    if ( ran[1] <= weig )
    {
      vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg;
      vflux[seg.slink[2].lx] = ffnew[0];
      vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ] = ffnew[1];

      nlink[ abs(seg.slink[2].lval) ]--;
      nlink[ abs(newvalseg) ]++;

      nflux[ff[0]]--;
      nflux[ff[1]]--;

      nflux[ffnew[0]]++;
      nflux[ffnew[1]]++;
        
      insert_next_segment = 0;
    }
  }
  else  // Insert next segment
  {
    seg.sdirmove = ran[1]<0.5; // 1=positive and 0=negative direction
    switch(seg.slink[2].lnu){
      case 0:
      rho++;
      break;

      case 1:
      rho += (rho!=0);
      break;

      case 2:
      rho = seg.slink[2].lnu - rho + (rho==0);
      break;
    }
    ranlxd(ran,2);
    seg.slink[3].lidx = int(ran[0]*3);
    int sseg = (seg.slink[2].lnu!=3)&(rho!=3);
    int xf[4];
    get_current_segments_and_links( rho, sseg, xf );

    // change sign of delta according to direction of move:
    int delta[5];
    if (rho==oldrho)
    {
      if (seg.sdirmove==oldsignseg) // go back
        a_eq_mb( delta, delta_seg );
      else
        a_eq_b( delta, delta_seg );
    }
    else
    {
      int var = (seg.slink[2].lnu<rho) + (seg.slink[2].lnu<oldrho);
      if (var==1)
      {
        if (seg.sdirmove!=oldsignseg)
          a_eq_mb( delta, delta_seg );
        else
          a_eq_b( delta, delta_seg );
      }
      else 
      {
        if (seg.sdirmove==oldsignseg)
          a_eq_mb( delta, delta_seg );
        else
          a_eq_b( delta, delta_seg );
      }
    }
    int newvalseg[3];
    newvalseg[0] = seg.slink[0].lval + delta[seg.slink[0].lidx];
    newvalseg[1] = seg.slink[1].lval + delta[seg.slink[1].lidx];
    newvalseg[2] = seg.svplaq + delta[4];

    int ll[2], ff[4], ffnew[4];
    ll[0] = vllink[seg.slink[0].lx][seg.slink[0].lnu];
    ll[1] = vllink[seg.slink[1].lx][seg.slink[1].lnu];

    if ( xf[3]<0 )
    {
      ff[0] = vflux[xf[0]];
      ff[1] = vflux[xf[1]];
      ff[2] = vflux[xf[2]];
      ff[3] = 0;

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval)
                       + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = 0;
    }
    else
    {
      ff[0] = vflux[xf[0]];
      ff[1] = vflux[xf[1]];
      ff[2] = vflux[xf[2]];
      ff[3] = vflux[xf[3]];

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = ff[3] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
    }
    // Check if fluxes > maximum
    ERROR_PLAQ(newvalseg[2]);
    ERROR_WORM(ffnew[0]); 
    ERROR_WORM(ffnew[1]);
    ERROR_WORM(ffnew[2]);
    ERROR_WORM(ffnew[3]);

    // calculate probability
    double weig = In[abs(newvalseg[2])]  - In[abs(seg.svplaq)]
               + Pn[ffnew[0]] - Pn[ff[0]]
               + Pn[ffnew[1]] - Pn[ff[1]]
               + Pn[ffnew[2]] - Pn[ff[2]]
               + Pn[ffnew[3]] - Pn[ff[3]]
               + fac[abs(seg.slink[0].lval)+ll[0]] - fac[abs(newvalseg[0])+ll[0]]
               + fac[abs(seg.slink[1].lval)+ll[1]] - fac[abs(newvalseg[1])+ll[1]];
    weig = exp( weig );

    if ( ran[1] <= weig )
    {
      vlink[seg.slink[0].lx][seg.slink[0].lnu] = newvalseg[0];
      vlink[seg.slink[1].lx][seg.slink[1].lnu] = newvalseg[1];
      vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg] = newvalseg[2];

      vflux[xf[0]] = ffnew[0];
      vflux[xf[1]] = ffnew[1];
      vflux[xf[2]] = ffnew[2];

      nlink[abs(seg.slink[0].lval)]--;
      nlink[abs(seg.slink[1].lval)]--;
      nplaq[abs(seg.svplaq)]--;

      nlink[abs(newvalseg[0])]++;
      nlink[abs(newvalseg[1])]++;
      nplaq[abs(newvalseg[2])]++;

      nflux[ff[0]]--;
      nflux[ff[1]]--;
      nflux[ff[2]]--;

      nflux[ffnew[0]]++;
      nflux[ffnew[1]]++;
      nflux[ffnew[2]]++;

      if ( xf[3]>=0 )
      {
        vflux[xf[3]] = ffnew[3];
        nflux[ff[3]]--;
        nflux[ffnew[3]]++;
      }

      oldsignseg = seg.ssignseg;
      if (seg.slink[3].lnu == rho){
        oldrho = seg.slink[2].lnu;
      }else{
        oldrho = rho;
      }
      // open_end --> open_init of next move
      seg.slink[2].lx = seg.slink[3].lx;
      seg.slink[2].lnu = seg.slink[3].lnu;
      seg.slink[2].lval = seg.slink[3].lval;
      oldlidx2 = seg.slink[3].lidx;

      a_eq_b( delta_seg, delta );
    }// accepted update
  }
}

//___________________________________________________________________________
void nworms(int nw)
{
  for(int iw=1; iw<=nw; iw++)
  {
    int inw = 0;
    while ( inw < NWORMS )
    {
      bool unsuccessful_worm = true;
      do{
        double ran;
        ranlxd( &ran, 1 );
        delta_seg[0] = 1-2*(ran<0.5);
        delta_seg[1] = delta_seg[0];
        delta_seg[2] = -delta_seg[1];
        delta_seg[3] = delta_seg[2];
        delta_seg[4] = delta_seg[3];
 
        insert_next_segment = 0;

        move0();
        while ( insert_next_segment == 1)
        {
          unsuccessful_worm = false;
          move1();
        };
        while ( insert_next_segment == 2 )
        {
          unsuccessful_worm = false;
          move();
        };
      }while(unsuccessful_worm);
      inw++;
    };

    // sweep of unscontrained variables
    sweep_l( );
  }
}

#endif
