#ifndef _SWEEPS_H
#define _SWEEPS_H

/**************************************************************************

	This file contains the following subroutins:

	void nsweeps(int ns)
	void metropolis_l()
	void metropolis_lbar()
	void metropolis_l_lbar()
	void metropolis_plaquette()
	void metropolis_s()
	void metropolis_sbar()
	void metropolis_s_sbar()
	void metropolis_meson()
	void metropolis_winding_loop()

**************************************************************************/

void metropolis_l(){

  int    ix,xnew,nu,lnew,nu2;
  int    delta,l,mnxx[2],mnxn;
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew = neib[ix][nu];
      nu2  = nu+nu;

      l  = dim[ix][nu2];
      delta = 3-6*(ran[nu]<0.5);
      lnew  = l + delta;
      
      if ( lnew >= 0 ){
  
      rho = vtau[delta+3]*facn[l]/facn[lnew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn    = mnxx[0] + delta;
      rho = rho*Tmn[mnxn][mnxx[1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[xnew][0];
      mnxx[1] = mnx[xnew][1];
      mnxn    = mnxx[1] + delta;
      rho = rho*Tmn[mnxx[0]][mnxn]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[nu+3] <= rho){
        dim[ix][nu2] = lnew;
        ndim        += delta;
        mnx[ix][0]  += delta;
        mnx[xnew][1] = mnxn;
      }
      }

    }
  }

}

//_________________________________________________________________________
void metropolis_lbar(){

  int    ix,xnew,nu,lbnew,nu2;
  int    delta,lb,mnxx[2],mnxn;
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew = neib[ix][nu];
      nu2  = nu+nu;

      lb = dim[ix][nu2+1];
      delta = 3-6*(ran[nu]<0.5);
      lbnew = lb + delta;
      
      if ( lbnew >= 0 ){
  
      rho = vtau[delta+3]*facn[lb]/facn[lbnew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn    = mnxx[1] + delta;
      rho = rho*Tmn[mnxx[0]][mnxn]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[xnew][0];
      mnxx[1] = mnx[xnew][1];
      mnxn    = mnxx[0] + delta;
      rho = rho*Tmn[mnxn][mnxx[1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[nu+3] <= rho){
        dim[ix][nu2+1] = lbnew;
        ndim        += delta;
        mnx[ix][1]  += delta;
        mnx[xnew][0] = mnxn;
      }
      }

    }
  }

}

//_________________________________________________________________________
void metropolis_l_lbar(){

  int    ix,xnew,nu,lnew,lbnew,nu2;
  int    delta,l,lb,mnxx[2],mnxn[2];
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew = neib[ix][nu];
      nu2  = nu+nu;

      l  = dim[ix][nu2];
      lb = dim[ix][nu2+1];
      delta = 1-2*(ran[nu]<0.5);
      lnew  = l + delta;
      lbnew = lb + delta;
      
      if ( (lnew>=0) && (lbnew>=0) ){
  
      rho = vtau[2*delta+3]*facn[l]/facn[lnew];
      rho = rho*facn[lb]/facn[lbnew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn[0] = mnxx[0] + delta;
      mnxn[1] = mnxx[1] + delta;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[xnew][0];
      mnxx[1] = mnx[xnew][1];
      mnxn[0] = mnxx[0] + delta;
      mnxn[1] = mnxx[1] + delta;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[nu+3] <= rho){
        dim[ix][nu2]   = lnew;
        dim[ix][nu2+1] = lbnew;
        ndim        += 2*delta;
        mnx[ix][0]  += delta;
        mnx[ix][1]  += delta;
        mnx[xnew][0] = mnxn[0];
        mnx[xnew][1] = mnxn[1]; 
      }
      }

    }
  }

}

//_________________________________________________________________________
void metropolis_plaquette(){

  int    ix,x[4],inu,nu0,nu1,l[2],lb[2],delta,lbnew[2],lnew[2];
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (inu=0; inu<3; inu++){
      nu0  = inu;
      nu1  = vnu[inu];

      x[0] = ix;
      x[1] = neib[ix][nu0];
      x[2] = neib[x[1]][nu1];
      x[3] = neib[ix][nu1];

      nu0 += nu0;
      nu1 += nu1;

      l[0]  = dim[x[0]][nu0];
      l[1]  = dim[x[0]][nu1];
      lb[0] = dim[x[1]][nu1+1];
      lb[1] = dim[x[3]][nu0+1];

      delta    = 1-2*(ran[inu]<0.5);
      lnew[0]  = l[0]  + delta;
      lnew[1]  = l[1]  - delta;
      lbnew[0] = lb[0] - delta;
      lbnew[1] = lb[1] + delta;

      if ( (lnew[0]>=0)  &&  (lnew[1]>=0)  &&  (lbnew[0]>=0)  &&  (lbnew[1]>=0) )
      {
        rho = facn[l[0]]/facn[lnew[0]];
        rho = rho*facn[l[1]]/facn[lnew[1]];
        rho = rho*facn[lb[0]]/facn[lbnew[0]];
        rho = rho*facn[lb[1]]/facn[lbnew[1]];

        if (ran[inu+3] <= rho)
        {
          dim[x[0]][nu0]   = lnew[0];
          dim[x[0]][nu1]   = lnew[1];
          dim[x[1]][nu1+1] = lbnew[0];
          dim[x[3]][nu0+1] = lbnew[1];
        }
      }
    }
  }

}

//_________________________________________________________________________
void metropolis_s(){

  int    ix,s,snew;
  int    delta,mnxx[2],mnxn;
  double ran[2*nsite],rho;

  ranlxd(ran,2*nsite);

  for (ix=0; ix<nsite; ix++)
  {
    s     = mon[ix][0];
    delta = 3-6*(ran[2*ix]<0.5);
    snew  = s + delta;
      
    if ( snew >= 0 ){

      rho = veta[delta+3]*facn[s]/facn[snew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn    = mnxx[0] + delta;
      rho = rho*Tmn[mnxn][mnxx[1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[2*ix+1] <= rho){
        mon[ix][0] = snew;
        nmon[0]   += delta;
        mnx[ix][0] = mnxn;
      }
    }

  }

}

//_________________________________________________________________________
void metropolis_sbar(){

  int    ix,sb,sbnew;
  int    delta,mnxx[2],mnxn;
  double ran[2*nsite],rho;

  ranlxd(ran,2*nsite);

  for (ix=0; ix<nsite; ix++)
  {
    sb    = mon[ix][1];
    delta = 3-6*(ran[2*ix]<0.5);
    sbnew = sb + delta;
      
    if ( sbnew >= 0 ){

      rho = vetabar[delta+3]*facn[sb]/facn[sbnew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn    = mnxx[1] + delta;
      rho = rho*Tmn[mnxx[0]][mnxn]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[2*ix+1] <= rho){
        mon[ix][1] = sbnew;
        nmon[1]   += delta;
        mnx[ix][1] = mnxn;
      }
    }

  }

}

//_________________________________________________________________________
void metropolis_s_sbar(){

  int    ix,s,sb,snew,sbnew;
  int    delta,mnxx[2],mnxn[2];
  double ran[2*nsite],rho;

  ranlxd(ran,2*nsite);

  for (ix=0; ix<nsite; ix++)
  {
    s     = mon[ix][0];
    sb    = mon[ix][1];
    delta = 1-2*(ran[2*ix]<0.5);
    snew  = s + delta;
    sbnew = sb + delta;
      
    if ( (sbnew>=0) && (snew>=0) ){

      rho = veta[delta+3]*facn[s]/facn[snew];
      rho = rho*vetabar[delta+3]*facn[sb]/facn[sbnew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn[0] = mnxx[0] + delta;
      mnxn[1] = mnxx[1] + delta;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[2*ix+1] <= rho){
        mon[ix][0] = snew;
        mon[ix][1] = sbnew;
        nmon[0]   += delta;
        nmon[1]   += delta;
        mnx[ix][0] = mnxn[0];
        mnx[ix][1] = mnxn[1];
      }
    }

  }

}


//_________________________________________________________________________
void metropolis_meson(){

  int    ix,xnew,nu,nu2,l,lnew;
  int    snew,sbnew,s,sb,delta;
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew  = neib[ix][nu];
      nu2   = nu+nu;
      delta = 1-2*(ran[nu]<0.5);

      l    = dim[ix][nu2];
      lnew = l + delta;

      s    = mon[ix][0];
      snew = s - delta;

      sb    = mon[xnew][1];
      sbnew = sb - delta;

      if ( (lnew>=0) && (snew>=0) && (sbnew>=0) )
      {
        rho = vtau[ delta+3]*facn[l]/facn[lnew];
        rho = rho*veta[-delta+3]*facn[s]/facn[snew];
        rho = rho*vetabar[-delta+3]*facn[sb]/facn[sbnew];

        if (ran[nu+3] <= rho)
        {
          dim[ix][nu2] = lnew;
          ndim        += delta;
          nmon[0]     -= delta;
          nmon[1]     -= delta;
          mon[ix][0]   = snew;
          mon[xnew][1] = sbnew;
        }
      }
    }
  }

}

//_________________________________________________________________________
void metropolis_winding_loop()
{
  int    i1,i2,i0,leng2,ix,xx,itype;
  int    delta,type,l,lnew,mnxx[2];
  double ran[3],rho;

  leng2 = leng*leng;

  for (itype=0; itype < 2; itype++){

  for ( i2=0; i2<leng; i2++ )
  {
    for ( i1=0; i1<leng; i1++ )
    {
       ranlxd(ran,2);
       ix    = i2*leng2 + i1*leng;
       delta = 1-2*(ran[0]<0.5);       
       type  = itype;//(ran[2]<0.5);
       rho   = 1.;

       for (i0=0; i0<leng; i0++)
       {
         xx   = ix + i0;
         l    = dim[xx][type];
         lnew = l + delta;
         if ( lnew < 0 ){
           rho = -1.;
           break;
         }
         else{
           rho     = rho*vtau[delta+3]*facn[l]/facn[lnew];
           mnxx[0] = mnx[xx][0];
           mnxx[1] = mnx[xx][1];
           rho = rho*Tmn[mnxx[0]+delta][mnxx[1]+delta]/Tmn[mnxx[0]][mnxx[1]];
         }
       }

       if ( ran[1] < rho ){
         for (i0=0; i0<leng; i0++){
           xx = ix + i0;
           dim[xx][type] += delta;
           ndim += delta;
           mnx[xx][0] += delta;
           mnx[xx][1] += delta;
         }
 
       }
    }
  }

  for ( i2=0; i2<leng; i2++ )
  {
    for ( i0=0; i0<leng; i0++ )
    {
       ix = i2*leng2 + i0;
       ranlxd(ran,2);
       delta = 1-2*(ran[0]<0.5);       
       type  = itype+2;//(ran[2]<0.5) + 2;
       rho   = 1.;
  
       for (i1=0; i1<leng; i1++){
         xx = ix + i1*leng;
         l  = dim[xx][type];
         lnew = l + delta;

         if ( lnew < 0 ){
           rho = -1.;
           break;
         }
         else{
           rho = rho*vtau[delta+3]*facn[l]/facn[lnew];
           mnxx[0] = mnx[xx][0];
           mnxx[1] = mnx[xx][1];
           rho = rho*Tmn[mnxx[0]+delta][mnxx[1]+delta]/Tmn[mnxx[0]][mnxx[1]];
         }
       }

       if ( ran[1] < rho ){
         for (i1=0; i1<leng; i1++)
         {
           xx = ix + i1*leng;
           dim[xx][type] += delta;
           ndim += delta;
           mnx[xx][0] += delta;
           mnx[xx][1] += delta;
         }
       }
    }
  }

  for ( i1=0; i1<leng; i1++ )
  {
    for ( i0=0; i0<leng; i0++ )
    {
       ix = i1*leng + i0;
       ranlxd(ran,2);
       delta = 1-2*(ran[0]<0.5);       
       type  = itype+4;//(ran[2]<0.5) + 4;
       rho   = 1.;

       for (i2=0; i2<leng; i2++){
         xx = ix + i2*leng2;
         l  = dim[xx][type];
         lnew = l + delta;
         if ( lnew < 0 ){
           rho = -1.;
           break;
         }
         else{
           rho = rho*vtau[delta+3]*facn[l]/facn[lnew];
           mnxx[0] = mnx[xx][0];
           mnxx[1] = mnx[xx][1];
           rho = rho*Tmn[mnxx[0]+delta][mnxx[1]+delta]/Tmn[mnxx[0]][mnxx[1]];
         }
       }

       if ( ran[1] < rho )
       {
         for (i2=0; i2<leng; i2++){
           xx = ix + i2*leng2;
           dim[xx][type] += delta;
           ndim += delta;
           mnx[xx][0] += delta;
           mnx[xx][1] += delta;
         }
       }
    }
  }

  }

}

//_________________________________________________________________________
void nsweeps(int ns)
{
  int is;
  
  for(is=1; is<=ns; is++)
  {
    metropolis_l();
    metropolis_lbar();
    metropolis_l_lbar();
    metropolis_plaquette();
#ifndef KAPPA0
    metropolis_s();
    metropolis_sbar();
    metropolis_s_sbar();
    metropolis_meson();
#endif
		metropolis_winding_loop();
  }
}

#endif
