#ifndef _WORMS_H
#define _WORMS_H

#ifdef SWA

//-----------------------------
// Data Structures
//-----------------------------

/* Link variable */
struct tlink
{
  int lidx; // position in the segment [0,3]
  int lval; // occupation number
  int lx;   // lattice site [0,nsites-1]
  int lnu;  // link direction [0,3]
};

/* Segment variable for worm updates */
struct segment
{
  tlink slink[4]; // 4 link variables: link[0] and link[1] = links to be modified 
                  //                   link[2] = current worm's head 
                  //                   link[3] = new worm's head
  int   svplaq;   // occupation number of plaquette variable
  int   sxplaq;   // site of the lower-left corner of the plaquette
  int   sdirmove; // if the segment will be inserted in positive or negative direction
  int   ssignseg; // if a positve or negative segment will be inserted
};

//-----------------------------
// Global Variables
//-----------------------------

int oldsignseg, oldrho, oldlidx2;
int delta_seg[5];
tlink link0;
segment seg;

//-----------------------------------------------------------------------------------------------------
//	Subroutines
//-----------------------------------------------------------------------------------------------------

#define ERROR(x,y) \
          if ( ((x)>=LENGTH) || ((y)>=LENGTH) ){ \
            cout << "FATAL ERROR: flux > max_LENGTH " << (x) << " " << (y) << endl; \
            return -1; \
          }

//-----------------------------------------------------------------------------------------------------
void get_current_segments_and_links( int rho, int sseg, int *xf, int ifield, int neib[][8] )
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
  seg.slink[0].lval = vlink[ seg.slink[0].lx ][ seg.slink[0].lnu + 4*ifield ];
  seg.slink[1].lval = vlink[ seg.slink[1].lx ][ seg.slink[1].lnu + 4*ifield ];
  seg.slink[3].lval = vlink[ seg.slink[3].lx ][ seg.slink[3].lnu + 4*ifield ];
  seg.svplaq = vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg];
}

//-----------------------------------------------------------------------------------------------------
inline void a_eq_b(int *a, int*b )
{
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
	a[3] = b[3];
	a[4] = b[4];
}

//-----------------------------------------------------------------------------------------------------
inline void a_eq_mb(int *a, int*b )
{
	a[0] = -b[0];
	a[1] = -b[1];
	a[2] = -b[2];
	a[3] = -b[3];
	a[4] = -b[4];
}

//-----------------------------------------------------------------------------------------------------
inline double seg_weight( int *newvalseg, int *ffnew, int *ff, int *ll )
{
    double weig = In[abs(newvalseg[2])]  - In[abs(seg.svplaq)]
               + Pn[ffnew[0]] - Pn[ff[0]]
               + Pn[ffnew[1]] - Pn[ff[1]]
               + Pn[ffnew[2]] - Pn[ff[2]]
               + Pn[ffnew[3]] - Pn[ff[3]]
               + fac[abs(seg.slink[0].lval)+ll[0]] - fac[abs(newvalseg[0])+ll[0]]
               + fac[abs(seg.slink[1].lval)+ll[1]] - fac[abs(newvalseg[1])+ll[1]];
    return exp( weig );
}

//-----------------------------------------------------------------------------------------------------
inline double link_weight( int *ffnew, int *ff, int ll, int newvalseg )
{
	double weig =   Pn[ffnew[0]] - Pn[ff[0]]
                + Pn[ffnew[1]] - Pn[ff[1]]
                + fac[ abs(seg.slink[2].lval) + ll ] - fac[ abs(newvalseg) + ll ];
  return exp( weig );
}

//-----------------------------------------------------------------------------------------------------
int move0( int ifield, int lnsite, int neib[][8] )
{
  double ran[4];
  ranlxd( ran, 4 );
  seg.slink[2].lx = int(ran[0]*lnsite);
  seg.slink[2].lnu = int(ran[1]*4);
  seg.slink[2].lval = vlink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ];

  // calculate probability
  link0.lval = 1-2*(ran[2]<0.5);
  int newvalseg = seg.slink[2].lval + link0.lval;

  int ll = vllink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ];
  int ff[2], ffnew[2];
  ff[0] = vflux[ seg.slink[2].lx ][ ifield ];
  ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ];

  ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
  ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);

	ERROR(ffnew[0],ffnew[1]);

  if ( ran[3] <= link_weight(ffnew,ff,ll, newvalseg) )
  {
    vlink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ] = newvalseg;
    vflux[ seg.slink[2].lx ][ ifield ] = ffnew[0];
    vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ] = ffnew[1];

    --nflux[ifield][ ff[0] ];
    --nflux[ifield][ ff[1] ];

    ++nflux[ifield][ ffnew[0] ];
    ++nflux[ifield][ ffnew[1] ];

    if ( seg.slink[2].lnu==3 )
    {
      nlink[ifield] += link0.lval;
    }

    return 1;
  }
  else
  {
    return 0;
  }
}

//-----------------------------------------------------------------------------------------------------
int move1( int ifield, int neib[][8] )
{
  double ran[2];
  ranlxd(ran,2);
  int rho = int(ran[0]*4); 

  if ( rho == 3 ) // insert link 
  {
    seg.slink[2].lval = vlink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ];
    int newvalseg = seg.slink[2].lval - link0.lval;

    int ll = vllink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ];
    int ff[2], ffnew[2];
    ff[0] = vflux[ seg.slink[2].lx ][ ifield ];
    ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ];

    ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
    ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);

		ERROR(ffnew[0],ffnew[1]);

    if ( ran[1] <= link_weight(ffnew,ff,ll,newvalseg) )
    {
      vlink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ] = newvalseg;
      vflux[ seg.slink[2].lx ][ ifield ] = ffnew[0];
      vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ] = ffnew[1];

      if ( seg.slink[2].lnu==3 )
      {
        nlink[ifield] -= link0.lval;
      }

      --nflux[ifield][ ff[0] ];
      --nflux[ifield][ ff[1] ];

      ++nflux[ifield][ ffnew[0] ];
      ++nflux[ifield][ ffnew[1] ];
        
      return 0;
    }
    else
    {
      return 1;
    }
  }
  else  // insert first segment
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
    get_current_segments_and_links( rho, sseg, xf, ifield, neib );

    // select sign:
    int delta[5];
    delta_seg[0] = -1;
    delta_seg[1] = -1;
    delta_seg[2] = +1;
    delta_seg[3] = +1;
    // ifield = 0 : link_0,link_1,link_2,link_3,pp = -1,-1,+1,+1,+1
    // ifield = 1 : link_0,link_1,link_2,link_3,pp = +1,+1,-1,-1,+1
    if ( ifield == 0 )
      delta_seg[4] = +1;
    else
      delta_seg[4] = -1;

    if ( link0.lval != delta_seg[seg.slink[2].lidx] )
			a_eq_mb( delta, delta_seg );
    else
			a_eq_b( delta, delta_seg );

    int newvalseg[3];
    newvalseg[0] = seg.slink[0].lval + delta[seg.slink[0].lidx];
    newvalseg[1] = seg.slink[1].lval + delta[seg.slink[1].lidx];
    newvalseg[2] = seg.svplaq + delta[4];

    int ll[2], ff[4], ffnew[4];
    ll[0] = vllink[ seg.slink[0].lx ][ seg.slink[0].lnu + 4*ifield ];
    ll[1] = vllink[ seg.slink[1].lx ][ seg.slink[1].lnu + 4*ifield ];

    if ( xf[3]<0 )
    {
      ff[0] = vflux[xf[0]][ifield];
      ff[1] = vflux[xf[1]][ifield];
      ff[2] = vflux[xf[2]][ifield];
      ff[3] = 0;

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval)
                       + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = 0;
    }
    else
    {
      ff[0] = vflux[xf[0]][ifield];
      ff[1] = vflux[xf[1]][ifield];
      ff[2] = vflux[xf[2]][ifield];
      ff[3] = vflux[xf[3]][ifield];

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = ff[3] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
    }
		ERROR(ffnew[0],ffnew[1]);
		ERROR(ffnew[2],ffnew[3]);

    // calculate probability
    if ( ran[1] <= seg_weight(newvalseg,ffnew,ff,ll) )
    {
      vlink[ seg.slink[0].lx ][seg.slink[0].lnu + 4*ifield ] = newvalseg[0];
      vlink[ seg.slink[1].lx ][seg.slink[1].lnu + 4*ifield ] = newvalseg[1];
      vplaq[ seg.sxplaq ][ seg.slink[2].lnu+rho-sseg ] = newvalseg[2];

      vflux[xf[0]][ifield] = ffnew[0];
      vflux[xf[1]][ifield] = ffnew[1];
      vflux[xf[2]][ifield] = ffnew[2];

      if ( seg.slink[0].lnu == 3 )
      {
        nlink[ifield] += delta[seg.slink[0].lidx]; 
      }
      if ( seg.slink[1].lnu == 3 )
      {
        nlink[ifield] += delta[seg.slink[1].lidx]; 
      }

      --nplaq[abs(seg.svplaq)];
      ++nplaq[abs(newvalseg[2])];

      --nflux[ifield][ff[0]];
      --nflux[ifield][ff[1]];
      --nflux[ifield][ff[2]];

      ++nflux[ifield][ffnew[0]];
      ++nflux[ifield][ffnew[1]];
      ++nflux[ifield][ffnew[2]];

      if ( xf[3]>=0 )
      {
        vflux[xf[3]][ifield] = ffnew[3];
        --nflux[ifield][ff[3]];
        ++nflux[ifield][ffnew[3]];
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

      return 2;
      // accepted update 
    }
    else
    {
      return 1;
    }
  }
}

//-----------------------------------------------------------------------------------------------------
int move( int ifield, int neib[][8] )
{
  double ran[2];
  ranlxd(ran,2);
  int rho = int(ran[0]*4); 

  if ( rho == 3 ) // try to close worm
  {
    int newvalseg = seg.slink[2].lval + delta_seg[oldlidx2];

    int ll = vllink[seg.slink[2].lx][seg.slink[2].lnu + 4*ifield];
    int ff[2], ffnew[2];
    ff[0] = vflux[ seg.slink[2].lx ][ ifield ];
    ff[1] = vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ];

    ffnew[0] = ff[0] + abs(newvalseg) - abs(seg.slink[2].lval);
    ffnew[1] = ff[1] + abs(newvalseg) - abs(seg.slink[2].lval);
    ERROR(ffnew[0],ffnew[1]);

    if ( ran[1] <= link_weight(ffnew,ff,ll,newvalseg) )
    {
      vlink[ seg.slink[2].lx ][ seg.slink[2].lnu + 4*ifield ] = newvalseg;
      vflux[ seg.slink[2].lx ][ ifield ] = ffnew[0];
      vflux[ neib[seg.slink[2].lx][seg.slink[2].lnu] ][ ifield ] = ffnew[1];

      if ( seg.slink[2].lnu == 3 )
      {
        nlink[ifield] += delta_seg[oldlidx2];
      }
      --nflux[ifield][ ff[0] ];
      --nflux[ifield][ ff[1] ];

      ++nflux[ifield][ ffnew[0] ];
      ++nflux[ifield][ ffnew[1]];
        
      return 0;
    }
    else
    {
      return 2;
    }
  }
  else  // Try to insert a new segment
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
    get_current_segments_and_links( rho, sseg, xf, ifield, neib );

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
    ll[0] = vllink[ seg.slink[0].lx ][ seg.slink[0].lnu + 4*ifield ];
    ll[1] = vllink[ seg.slink[1].lx ][ seg.slink[1].lnu + 4*ifield ];

    if ( xf[3]<0 )
    {
      ff[0] = vflux[xf[0]][ifield];
      ff[1] = vflux[xf[1]][ifield];
      ff[2] = vflux[xf[2]][ifield];
      ff[3] = 0;

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval)
                       + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = 0;
    }
    else
    {
      ff[0] = vflux[xf[0]][ifield];
      ff[1] = vflux[xf[1]][ifield];
      ff[2] = vflux[xf[2]][ifield];
      ff[3] = vflux[xf[3]][ifield];

      ffnew[0] = ff[0] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[1] = ff[1] + abs(newvalseg[0]) - abs(seg.slink[0].lval);
      ffnew[2] = ff[2] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
      ffnew[3] = ff[3] + abs(newvalseg[1]) - abs(seg.slink[1].lval);
    }

   	ERROR(ffnew[0],ffnew[1]);
		ERROR(ffnew[2],ffnew[3]);

    // calculate probability
    if ( ran[1] <= seg_weight(newvalseg,ffnew,ff,ll) )
    {
      vlink[ seg.slink[0].lx ][ seg.slink[0].lnu + 4*ifield ] = newvalseg[0];
      vlink[ seg.slink[1].lx ][ seg.slink[1].lnu + 4*ifield ] = newvalseg[1];
      vplaq[ seg.sxplaq ][ seg.slink[2].lnu+rho-sseg ] = newvalseg[2];

      vflux[xf[0]][ifield] = ffnew[0];
      vflux[xf[1]][ifield] = ffnew[1];
      vflux[xf[2]][ifield] = ffnew[2];

      if ( seg.slink[0].lnu == 3 ) 
      {
        nlink[ifield] += delta[seg.slink[0].lidx];
      }
      if ( seg.slink[1].lnu == 3 )
      {
        nlink[ifield] += delta[seg.slink[1].lidx];
      }

      --nplaq[abs(seg.svplaq)];
      ++nplaq[abs(newvalseg[2])];

      --nflux[ifield][ ff[0] ];
      --nflux[ifield][ ff[1] ];
      --nflux[ifield][ ff[2] ];

      ++nflux[ifield][ ffnew[0] ];
      ++nflux[ifield][ ffnew[1] ];
      ++nflux[ifield][ ffnew[2] ];

      if ( xf[3]>=0 )
      {
        vflux[xf[3]][ifield] = ffnew[3];
        --nflux[ifield][ff[3] ];
        ++nflux[ifield][ffnew[3] ];
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
  
    return 2;
  }
}

//-----------------------------------------------------------------------------------------------------
int surface_worms( int ifield, int lnsite, int neib[][8] )
{
  int next_segment = move0( ifield, lnsite, neib );

  if ( next_segment == -1 )
  {
    return -1;
  }
  else if ( next_segment == 0 )
  { 
    return 0;
  }
  else
  {
    while ( next_segment == 1 ) next_segment = move1( ifield, neib );
    while ( next_segment == 2 ) next_segment = move( ifield, neib );

     return 1;
  }

}

#endif

#endif 
