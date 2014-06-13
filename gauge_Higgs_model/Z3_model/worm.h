#ifndef _WORM_H
#define _WORM_H

/************************************************************************* 
*		This file contains:
*
*		Total update:
*		void nworms(int nw)
*
*		Worm update for constrained variables:
*		void move0()
*		void move1()
*		void move()
*		void get_current_segments_and_links( int rho, int sseg, int *xf )
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

// STRUCT VARS
struct tlink
{
  /* Link variable */
  int lidx;// position in the segment [0,3]
  int lval; // occupation number: +1,0,-1
  int lx;   // lattice site [0,nsites-1]
  int lnu;  // link direction [0,3]
};

struct segment
{
  /* Segment variable for worm updates */
  tlink slink[4]; // 4 link variables: link[0] and link[1] = links to be modified 
                  //                   link[2] = current worm's head 
                  //                   link[3] = new worm's head
  int svplaq;   // occupation number of plaquette variable
  int sxplaq;   // site of the lower-left corner of the plaquette
  int sdirmove; // if the segment will be inserted in positive or negative direction
  int ssignseg; // if a positve or negative segment will be inserted
};

// Global variables

//  triadd[0][0] = 2;
//  triadd[0][1] = 0;
//  triadd[0][2] = 1;
//  triadd[1][0] = 0;
//  triadd[1][1] = 1;
//  triadd[1][2] = 2;
//  triadd[2][0] = 1;  
//  triadd[2][1] = 2;
//  triadd[2][2] = 0;
int triadd[3][3] = {{2,0,1},{0,1,2},{1,2,0}};
int insert_next_segment;   // if 0 then worm terminates
int oldsignseg,oldrho,real_world,oldlidx2;   // keep track of old values 
int delta_seg[5] = {0,0,2,2,2}; //  link_0,link_1,link_2,link_3,pp = -1,-1,+1,+1,+1
tlink  link0;   // link where the worm starts
segment seg;

  
// Subroutines
//________________________________________________________________________
void get_current_segments_and_links( int rho, int sseg )
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
        break;        
      }      
    }
  }
  seg.slink[0].lval = vlink[seg.slink[0].lx][seg.slink[0].lnu];
  seg.slink[1].lval = vlink[seg.slink[1].lx][seg.slink[1].lnu];
  seg.slink[3].lval = vlink[seg.slink[3].lx][seg.slink[3].lnu];
  seg.svplaq = vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg];

}

//________________________________________________________________________
inline void a_eq_b(int *a, int *b )
{
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
	a[3] = b[3];
	a[4] = b[4];
}

//________________________________________________________________________
inline void a_eq_mb(int *a, int *b )
{
	a[0] = 2-b[0];
	a[1] = 2-b[1];
	a[2] = 2-b[2];
	a[3] = 2-b[3];
	a[4] = 2-b[4];
}

//_________________________________________________________________________
void move0()
{
	//
	// Insert first link
	//
  double weig,ran[4];

  ranlxd(ran,4);
  seg.slink[2].lx  = int(ran[0]*nsite);
  seg.slink[2].lnu = int(ran[1]*4);
  seg.slink[2].lval = vlink[seg.slink[2].lx][seg.slink[2].lnu];

  /* choose value of first link randomly */
  link0.lval = 2*(ran[2]<0.5);
  int newvalseg  = triadd[ seg.slink[2].lval ][ link0.lval ];

  /* calculate Metropolis probability */
  if ( seg.slink[2].lnu!=3 ){
    weig = logbb*( abs(newvalseg-1) - abs(seg.slink[2].lval-1) );
  }else{
    weig = link4weight[newvalseg] - link4weight[seg.slink[2].lval];
  }

  /* accept/reject step */
  if ( ran[3]<=exp(weig) )
  {
    insert_next_segment = 1;
    
    vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg;
    nlink[seg.slink[2].lnu][seg.slink[2].lval]--;
    nlink[seg.slink[2].lnu][newvalseg]++;
  }
}

//_________________________________________________________________________
void move1()
{
	// 
	// The worm tries to insert the first segment
	//
  int sseg,newvalseg[3],delta[5];
  double weig,ran[2];

  ranlxd(ran,2);
  int rho = int(ran[0]*4); 

  if ( rho == 3 )// try to close the worm undoing move0()
  {
    seg.slink[2].lval = vlink[seg.slink[2].lx][seg.slink[2].lnu];
    newvalseg[0] = triadd[ seg.slink[2].lval ][ 2-link0.lval ];

    /* Metropolis probability */
    if ( seg.slink[2].lnu!=3 ){
      weig = logbb*( abs(newvalseg[0]-1) - abs(seg.slink[2].lval-1) );
    }else{
      weig = link4weight[newvalseg[0]] - link4weight[seg.slink[2].lval];
    }
    /* accept/reject step */
    if ( ran[1]<=exp(weig) )
    {
      vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg[0];
      nlink[seg.slink[2].lnu][ seg.slink[2].lval ]--;
      nlink[seg.slink[2].lnu][ newvalseg[0] ]++;
        
      insert_next_segment = 0;
    }

  }
  else // Insert segment
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
    sseg = (seg.slink[2].lnu!=3)&(rho!=3);
    get_current_segments_and_links( rho, sseg );

    // select sign:
    if ( link0.lval != delta_seg[seg.slink[2].lidx] )
			a_eq_mb( delta, delta_seg );
    else
			a_eq_b( delta, delta_seg );

    newvalseg[0] = triadd[ seg.slink[0].lval ][ delta[seg.slink[0].lidx] ];
    newvalseg[1] = triadd[ seg.slink[1].lval ][ delta[seg.slink[1].lidx] ];
    newvalseg[2] = triadd[ seg.svplaq ][ delta[4] ];

    // calculate probability
    weig = logbbt*( abs(newvalseg[2]-1) - abs(seg.svplaq-1) );

    if ( seg.slink[0].lnu!=3 ){
     weig += logbb*( abs(newvalseg[0]-1) - abs(seg.slink[0].lval-1) );
    }else{
      weig += link4weight[newvalseg[0]] - link4weight[seg.slink[0].lval];
    }
    if ( seg.slink[1].lnu!=3 ){
      weig += logbb*( abs(newvalseg[1]-1) - abs(seg.slink[1].lval-1) );
    }else{
      weig += link4weight[newvalseg[1]] - link4weight[seg.slink[1].lval];
    }

    /* accept/reject step */
    if (ran[1]<=exp(weig))
    {
      insert_next_segment = 2;

      vlink[seg.slink[0].lx][seg.slink[0].lnu] = newvalseg[0];
      vlink[seg.slink[1].lx][seg.slink[1].lnu] = newvalseg[1];
      vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg] = newvalseg[2];

      nlink[seg.slink[0].lnu][seg.slink[0].lval]--;
      nlink[seg.slink[1].lnu][seg.slink[1].lval]--;
      nplaq[seg.svplaq]--;

      nlink[seg.slink[0].lnu][newvalseg[0]]++;
      nlink[seg.slink[1].lnu][newvalseg[1]]++;
      nplaq[newvalseg[2]]++;    

      oldsignseg = seg.ssignseg;
      if (seg.slink[3].lnu == rho){
        oldrho = seg.slink[2].lnu;
      }else{
        oldrho = rho;
      }
      /* new position of worm's head --> current worm's head position */
      seg.slink[2].lx = seg.slink[3].lx;
      seg.slink[2].lnu = seg.slink[3].lnu;
      seg.slink[2].lval = seg.slink[3].lval;
      oldlidx2 = seg.slink[3].lidx;

			a_eq_b( delta_seg, delta );

    }// accepted update
  }
}

//_________________________________________________________________________
void move()
{
  int    sseg,rho,newvalseg[3],delta[5],var;
  double weig,ran[2];

  ranlxd(ran,2);
  rho = int(ran[0]*4); 

  if ( rho == 3 ) // try to close worm
  {
    newvalseg[0] = triadd[ seg.slink[2].lval ][ delta_seg[oldlidx2] ];

    /* calculate Metropolis probability */
    if ( seg.slink[2].lnu!=3 ){
      weig = logbb*( abs(newvalseg[0]-1) - abs(seg.slink[2].lval-1) );
    }else{
      weig = link4weight[newvalseg[0]] - link4weight[seg.slink[2].lval];
    }
    /* accept/reject step */
    if ( ran[1]<=exp(weig) )
    {
      vlink[seg.slink[2].lx][seg.slink[2].lnu] = newvalseg[0];
      nlink[seg.slink[2].lnu][ seg.slink[2].lval ]--;
      nlink[seg.slink[2].lnu][ newvalseg[0] ]++;
        
      insert_next_segment = false;
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
    sseg = (seg.slink[2].lnu!=3)&(rho!=3);
    get_current_segments_and_links( rho, sseg );

    // change sign of delta according to direction of move:
    if (rho==oldrho)
    {
      if (seg.sdirmove==oldsignseg) /* worm visits last updated segment */
				a_eq_mb( delta, delta_seg ); 
      else
				a_eq_b( delta, delta_seg );
    }
    else
    {
      var = (seg.slink[2].lnu<rho) + (seg.slink[2].lnu<oldrho);
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
    newvalseg[0] = triadd[ seg.slink[0].lval ][ delta[seg.slink[0].lidx] ];
    newvalseg[1] = triadd[ seg.slink[1].lval ][ delta[seg.slink[1].lidx] ];
    newvalseg[2] = triadd[ seg.svplaq ][ delta[4] ];

    // calculate probability
    weig = logbbt*( abs(newvalseg[2]-1) - abs(seg.svplaq-1) );

    if ( seg.slink[0].lnu!=3 ){
     weig += logbb*( abs(newvalseg[0]-1) - abs(seg.slink[0].lval-1) );
    }else{
      weig += link4weight[newvalseg[0]] - link4weight[seg.slink[0].lval];
    }
    if ( seg.slink[1].lnu!=3 ){
      weig += logbb*( abs(newvalseg[1]-1) - abs(seg.slink[1].lval-1) );
    }else{
      weig += link4weight[newvalseg[1]] - link4weight[seg.slink[1].lval];
    }

    /* accept/reject step */
    if (ran[1]<=exp(weig))
    {
      vlink[seg.slink[0].lx][seg.slink[0].lnu] = newvalseg[0];
      vlink[seg.slink[1].lx][seg.slink[1].lnu] = newvalseg[1];
      vplaq[seg.sxplaq][seg.slink[2].lnu+rho-sseg] = newvalseg[2];

      nlink[seg.slink[0].lnu][seg.slink[0].lval]--;
      nlink[seg.slink[1].lnu][seg.slink[1].lval]--;
      nplaq[seg.svplaq]--;

      nlink[seg.slink[0].lnu][newvalseg[0]]++;
      nlink[seg.slink[1].lnu][newvalseg[1]]++;
      nplaq[newvalseg[2]]++;    

      oldsignseg = seg.ssignseg;
      if (seg.slink[3].lnu == rho){
        oldrho = seg.slink[2].lnu;
      }else{
        oldrho = rho;
      }
      /* new position of worm's head --> current worm's head position */
      seg.slink[2].lx = seg.slink[3].lx;
      seg.slink[2].lnu = seg.slink[3].lnu;
      seg.slink[2].lval = seg.slink[3].lval;
      oldlidx2 = seg.slink[3].lidx;

			a_eq_b( delta_seg, delta );
    }// accepted update
  }
}

//________________________________________________________________________
void nworms(int nw)
{
  bool unsuccessful_worm;
  double ran[1];

  for(int iw=1; iw<=nw; iw++)
  {
    int inw= 0;
    while(inw < nsite)
    {
	    unsuccessful_worm = true;
  	  do{
  	    ranlxd(ran,1);
  	    delta_seg[0] = 2*(ran[0]<0.5);
  	    delta_seg[1] = delta_seg[0];
  	    delta_seg[2] = 2-delta_seg[1];
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
  }
}

#endif
