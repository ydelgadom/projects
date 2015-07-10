#ifndef _SWEEPS_H
#define _SWEEPS_H

/*********************************************************************

  Subroutines in this file:

  void metropolis_k()
  void metropolis_kbar()
  void metropolis_plaquette()
  void metropolis_s()
  void metropolis_sbar()
  void metropolis_meson()
  void metropolis_winding_loop()

**********************************************************************/
void metropolis_k(){

  int    ix,xnew,nu,knew,abskb,nu2;
  int    delta,k,mnxx[2],mnxn[2];
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew = neib[ix][nu];
      nu2  = nu+nu;

      abskb = int(abs(float(dim[ix][nu2])));
      k     = dim[ix][nu2+1];
      delta = 1-2*(ran[nu]<0.5);
      knew  = k + delta;
      
      if ( knew >= 0 ){
  
      rho = vtau[2*delta+7]*facn[abskb + k]/facn[abskb + knew];
      rho = rho*facn[k]/facn[knew];

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
        dim[ix][nu2+1] = knew;
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
void metropolis_kbar(){

  int    ix,xnew,nu,nu2,kbnew,abskb,abskbnew;
  int    delta,k,diff,mnxx[2],mnxn[2];
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew = neib[ix][nu];
      nu2  = nu+nu;

      diff  = dim[ix][nu2];
      k     = dim[ix][nu2+1];
      abskb = int(abs(float(diff)));
      delta = 3-6*(ran[nu]<0.5);
      kbnew = diff + delta;
      abskbnew = int(abs(float(kbnew)));
      diff  = abskbnew - abskb;
 
      rho = vtau[diff+7]*facn[abskb + k]/facn[abskbnew + k];

      abskbnew = (diff + delta)/2;
      abskb    = (diff - delta)/2;

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn[0] = mnxx[0] + abskbnew;
      mnxn[1] = mnxx[1] + abskb;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[xnew][0];
      mnxx[1] = mnx[xnew][1];
      mnxn[0] = mnxx[0] + abskb;
      mnxn[1] = mnxx[1] + abskbnew;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]]; 

      if (ran[nu+3] <= rho){
        dim[ix][nu2] = kbnew;
        ndim        += diff;
        mnx[ix][0]  += abskbnew;
        mnx[ix][1]  += abskb;
        mnx[xnew][0] = mnxn[0];
        mnx[xnew][1] = mnxn[1];
      }

    }
  }

}

//_________________________________________________________________________
void metropolis_plaquette(){

  int    ix,x[4],inu,nu0,nu1,abskb[4],abskbnew[4],k[4];
  int    delta,mnxn[4][2],kb[4],kbnew[4],mnxx[2];
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

      k[0] = dim[x[0]][nu0+1];
      k[1] = dim[x[1]][nu1+1];
      k[2] = dim[x[3]][nu0+1];
      k[3] = dim[x[0]][nu1+1];

      kb[0] = dim[x[0]][nu0];
      kb[1] = dim[x[1]][nu1];
      kb[2] = dim[x[3]][nu0];
      kb[3] = dim[x[0]][nu1];

      abskb[0] = int(abs(float(kb[0])));
      abskb[1] = int(abs(float(kb[1])));
      abskb[2] = int(abs(float(kb[2])));
      abskb[3] = int(abs(float(kb[3])));

      delta    = 1-2*(ran[inu]<0.5);
      kbnew[0] = kb[0] + delta;
      kbnew[1] = kb[1] + delta;
      kbnew[2] = kb[2] - delta;
      kbnew[3] = kb[3] - delta;

      abskbnew[0] = int(abs(float(kbnew[0])));
      abskbnew[1] = int(abs(float(kbnew[1])));
      abskbnew[2] = int(abs(float(kbnew[2])));
      abskbnew[3] = int(abs(float(kbnew[3])));

      kb[0] = abskbnew[0] - abskb[0];
      kb[1] = abskbnew[1] - abskb[1];
      kb[2] = abskbnew[2] - abskb[2];
      kb[3] = abskbnew[3] - abskb[3];

      delta = kb[0] + kb[1] + kb[2] + kb[3];

      rho = vtau[delta+7];
      rho = rho*facn[abskb[0] + k[0]]/facn[abskbnew[0] + k[0]];
      rho = rho*facn[abskb[1] + k[1]]/facn[abskbnew[1] + k[1]];
      rho = rho*facn[abskb[2] + k[2]]/facn[abskbnew[2] + k[2]];
      rho = rho*facn[abskb[3] + k[3]]/facn[abskbnew[3] + k[3]];

      mnxx[0] = mnx[x[0]][0];
      mnxx[1] = mnx[x[0]][1];
      mnxn[0][0] = mnxx[0] + (kb[0] + kb[3])/2;
      mnxn[0][1] = mnxx[1] + (kb[0] + kb[3])/2;
      rho = rho*Tmn[mnxn[0][0]][mnxn[0][1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[x[1]][0];
      mnxx[1] = mnx[x[1]][1];
      mnxn[1][0] = mnxx[0] + (kb[0] + kb[1])/2;
      mnxn[1][1] = mnxx[1] + (kb[0] + kb[1])/2;
      rho = rho*Tmn[mnxn[1][0]][mnxn[1][1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[x[2]][0];
      mnxx[1] = mnx[x[2]][1];
      mnxn[2][0] = mnxx[0] + (kb[1] + kb[2])/2;
      mnxn[2][1] = mnxx[1] + (kb[1] + kb[2])/2;
      rho = rho*Tmn[mnxn[2][0]][mnxn[2][1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[x[3]][0];
      mnxx[1] = mnx[x[3]][1];
      mnxn[3][0] = mnxx[0] + (kb[2] + kb[3])/2;
      mnxn[3][1] = mnxx[1] + (kb[2] + kb[3])/2;
      rho = rho*Tmn[mnxn[3][0]][mnxn[3][1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[inu+3] <= rho){
        dim[x[0]][nu0] = kbnew[0];
        dim[x[1]][nu1] = kbnew[1];
        dim[x[3]][nu0] = kbnew[2];
        dim[x[0]][nu1] = kbnew[3];

        ndim += delta;

        mnx[x[0]][0] = mnxn[0][0];
        mnx[x[0]][1] = mnxn[0][1];
        mnx[x[1]][0] = mnxn[1][0];
        mnx[x[1]][1] = mnxn[1][1];
        mnx[x[2]][0] = mnxn[2][0];
        mnx[x[2]][1] = mnxn[2][1];
        mnx[x[3]][0] = mnxn[3][0];
        mnx[x[3]][1] = mnxn[3][1];
      }

    }
  }

}
//_________________________________________________________________________
void metropolis_s(){

  int    ix,snew,abssb;
  int    delta,s,mnxx[2],mnxn[2];
  double ran[2*nsite],rho;

  ranlxd(ran,2*nsite);

  for (ix=0; ix<nsite; ix++)
  {
    s     = mon[ix][0];
    abssb = int(abs(float(mon[ix][1])));
    delta = 1-2*(ran[2*ix]<0.5);
    snew  = s + delta;
      
    if ( snew >= 0 ){

      rho = vkappa[2*delta+7]*facn[abssb + s]/facn[abssb + snew];
      rho = rho*facn[s]/facn[snew];

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn[0] = mnxx[0] + delta;
      mnxn[1] = mnxx[1] + delta;
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      if (ran[2*ix+1] <= rho){
        mon[ix][0] = snew;
        nmon[0]   += 2*delta;
        mnx[ix][0] = mnxn[0];
        mnx[ix][1] = mnxn[1];
      }
    }

  }

}

//_________________________________________________________________________
void metropolis_sbar(){

  int    ix,sbnew,abssb,abssbnew;
  int    delta,s,sb,mnxx[2],mnxn[2];
  double ran[2],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,2);

    s   = mon[ix][0];
    sb  = mon[ix][1];
    abssb = int(abs(float(sb)));
    delta = 3-6*(ran[0]<0.5);
    sbnew = sb + delta;
    abssbnew = int(abs(float(sbnew)));
    sb  = abssbnew - abssb; 
 
    rho = vkappa[sb+7]*vemu[delta+7]*facn[abssb + s]/facn[abssbnew + s];

    abssbnew = (sb + delta)/2;
    abssb    = (sb - delta)/2;

    mnxx[0] = mnx[ix][0];
    mnxx[1] = mnx[ix][1];
    mnxn[0] = mnxx[0] + abssbnew;
    mnxn[1] = mnxx[1] + abssb;
    rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

    if (ran[1] <= rho){
      mon[ix][1] = sbnew;
      nmon[1]   += delta;
      nmon[2]   += sb;
      mnx[ix][0] = mnxn[0];
      mnx[ix][1] = mnxn[1];
    }

  }

}

//_________________________________________________________________________
void metropolis_meson(){

  int    ix,xnew,nu,nu2,kbnew,abskb,abskbnew;
  int    sbnew[2],abssb[2],abssbnew[2],s[2],sb[2];
  int    delta,k,kb,mnxx[2],mnxn[2];
  double ran[6],rho;

  for (ix=0; ix<nsite; ix++)
  {
    ranlxd(ran,6);

    for (nu=0; nu<3; nu++){
      xnew  = neib[ix][nu];
      nu2   = nu+nu;
      delta = 1-2*(ran[nu]<0.5);

      kb    = dim[ix][nu2];
      k     = dim[ix][nu2+1];
      abskb = int(abs(float(kb)));
      kbnew = kb + delta;
      abskbnew = int(abs(float(kbnew)));
      kb    = abskbnew - abskb;
 
      rho = vtau[kb+7]*facn[abskb + k]/facn[abskbnew + k];

      abskbnew = (kb + delta)/2;
      abskb    = (kb - delta)/2;

      s[0]   = mon[ix][0];
      sb[0]  = mon[ix][1];
      abssb[0] = int(abs(float(sb[0])));
      sbnew[0] = sb[0] - delta;
      abssbnew[0] = int(abs(float(sbnew[0])));
      sb[0]  = abssbnew[0] - abssb[0]; 

      s[1]   = mon[xnew][0];
      sb[1]  = mon[xnew][1];
      abssb[1] = int(abs(float(sb[1])));
      sbnew[1] = sb[1] + delta;
      abssbnew[1] = int(abs(float(sbnew[1])));
      sb[1]  = abssbnew[1] - abssb[1];
 
      rho = rho*vkappa[sb[0]+sb[1]+7]*facn[abssb[0] + s[0]]/facn[abssbnew[0] + s[0]];
      rho = rho*facn[abssb[1] + s[1]]/facn[abssbnew[1] + s[1]];

      abssbnew[0] = (sb[0] - delta)/2;
      abssb[0]    = (sb[0] + delta)/2;

      abssbnew[1] = (sb[1] + delta)/2;
      abssb[1]    = (sb[1] - delta)/2;

      mnxx[0] = mnx[ix][0];
      mnxx[1] = mnx[ix][1];
      mnxn[0] = mnxx[0] + abskbnew + abssbnew[0];
      mnxn[1] = mnxx[1] + abskb + abssb[0];
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

      mnxx[0] = mnx[xnew][0];
      mnxx[1] = mnx[xnew][1];
      mnxn[0] = mnxx[0] + abskb + abssbnew[1];
      mnxn[1] = mnxx[1] + abskbnew + abssb[1];
      rho = rho*Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]]; 

      if (ran[nu+3] <= rho){
        dim[ix][nu2] = kbnew;
        ndim        += kb;
        mon[ix][1]   = sbnew[0];
        mon[xnew][1] = sbnew[1];
        nmon[2]     += (sb[0] + sb[1]);

        mnx[ix][0]  += abskbnew + abssbnew[0];
        mnx[ix][1]  += abskb + abssb[0];
        mnx[xnew][0] = mnxn[0];
        mnx[xnew][1] = mnxn[1];
      }

    }
  }

}

//_________________________________________________________________________
void metropolis_winding_loop()
{

  int    i1,i2,i3,leng2,ix,xx,xnew;
  int    delta,diff,mnxx[2],mnxn[2];
  int    k,kb,kbnew,abskb,abskbnew;
  double ran[2],rho;

  leng2 = leng*leng;

  // loop in z direction
  for (i1=0; i1<leng; i1++)
  {
    for (i2=0; i2<leng; i2++)
    {
      ranlxd(ran,2);
      delta = 1-2*(ran[0]<0.5);

      ix   = i2*leng + i1;
      diff = 0;
      rho  = 1.;
      for (i3=0; i3<leng; i3++)
      {
        xx   = ix + i3*leng2;
        xnew = xx + leng2;

        kb = dim[xx][4];
        k  = dim[xx][5];
        kbnew = kb + delta;
        abskb = abs(kb);
        abskbnew = abs(kbnew);
        rho *= facn[abskb + k]/facn[abskbnew + k];
        k = abskbnew - abskb;
        diff += k;
        rho *= vtau[k+7];

        abskbnew = (k + delta)/2;
        abskb    = (k - delta)/2;

        mnxx[0] = mnx[xx][0];
        mnxx[1] = mnx[xx][1];
        mnxn[0] = mnxx[0] + abskbnew;
        mnxn[1] = mnxx[1] + abskb;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

        mnxx[0] = mnx[xnew][0];
        mnxx[1] = mnx[xnew][1];
        mnxn[0] = mnxx[0] + abskb;
        mnxn[1] = mnxx[1] + abskbnew;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]]; 
      }

      if (ran[1] <= rho)
      {
        ndim += diff;
        for (i3=0; i3<leng; i3++)
        {
          xx   = ix + i3*leng2;
          xnew = xx + leng2;

          k = dim[xx][4];
          dim[xx][4] = k + delta;
          k = abs(k+delta) - abs(k);
          abskbnew = (k + delta)/2;
          abskb    = (k - delta)/2;
          mnx[xx][0] += abskbnew;
          mnx[xx][1] += abskb;
          mnx[xnew][0] += abskb;
          mnx[xnew][1] += abskbnew;
        }
      }
    }
  }

  // loop in x direction
  for (i1=0; i1<leng; i1++)
  {
    for (i2=0; i2<leng; i2++)
    {
      ranlxd(ran,2);
      delta = 1-2*(ran[0]<0.5);

      ix = i2*leng2 + i1*leng;
      diff = 0;
      rho  = 1.;
      for (i3=0; i3<leng; i3++)
      {
        xx   = ix + i3;
        xnew = xx + 1;

        kb = dim[xx][0];
        k  = dim[xx][1];
        kbnew = kb + delta;
        abskb = abs(kb);
        abskbnew = abs(kbnew);
        rho *= facn[abskb + k]/facn[abskbnew + k];
        k = abskbnew - abskb;
        diff += k;
        rho *= vtau[k+7];

        abskbnew = (k + delta)/2;
        abskb    = (k - delta)/2;

        mnxx[0] = mnx[xx][0];
        mnxx[1] = mnx[xx][1];
        mnxn[0] = mnxx[0] + abskbnew;
        mnxn[1] = mnxx[1] + abskb;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

        mnxx[0] = mnx[xnew][0];
        mnxx[1] = mnx[xnew][1];
        mnxn[0] = mnxx[0] + abskb;
        mnxn[1] = mnxx[1] + abskbnew;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]]; 
      }

      if (ran[1] <= rho)
      {
        ndim += diff;
        for (i3=0; i3<leng; i3++)
        {
          xx   = ix + i3;
          xnew = xx + 1;

          k = dim[xx][0];
          dim[xx][0] = k + delta;
          k = abs(k+delta) - abs(k);
          abskbnew = (k + delta)/2;
          abskb    = (k - delta)/2;
          mnx[xx][0] += abskbnew;
          mnx[xx][1] += abskb;
          mnx[xnew][0] += abskb;
          mnx[xnew][1] += abskbnew;
        }
      }
    }
  }

    // loop in y direction
  for (i1=0; i1<leng; i1++)
  {
    for (i2=0; i2<leng; i2++)
    {
      ranlxd(ran,2);
      delta = 1-2*(ran[0]<0.5);

      ix   = i2*leng2 + i1;
      diff = 0;
      rho  = 1.;
      for (i3=0; i3<leng; i3++)
      {
        xx   = ix + i3*leng;
        xnew = xx + leng;

        kb = dim[xx][2];
        k  = dim[xx][3];
        kbnew = kb + delta;
        abskb = abs(kb);
        abskbnew = abs(kbnew);
        rho *= facn[abskb + k]/facn[abskbnew + k];
        k = abskbnew - abskb;
        diff += k;
        rho *= vtau[k+7];

        abskbnew = (k + delta)/2;
        abskb    = (k - delta)/2;

        mnxx[0] = mnx[xx][0];
        mnxx[1] = mnx[xx][1];
        mnxn[0] = mnxx[0] + abskbnew;
        mnxn[1] = mnxx[1] + abskb;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]];

        mnxx[0] = mnx[xnew][0];
        mnxx[1] = mnx[xnew][1];
        mnxn[0] = mnxx[0] + abskb;
        mnxn[1] = mnxx[1] + abskbnew;
        rho *= Tmn[mnxn[0]][mnxn[1]]/Tmn[mnxx[0]][mnxx[1]]; 
      }

      if (ran[1] <= rho)
      {
        ndim += diff;
        for (i3=0; i3<leng; i3++)
        {
          xx   = ix + i3*leng;
          xnew = xx + leng;

          k = dim[xx][2];
          dim[xx][2] = k + delta;
          k = abs(k+delta) - abs(k);
          abskbnew = (k + delta)/2;
          abskb    = (k - delta)/2;
          mnx[xx][0] += abskbnew;
          mnx[xx][1] += abskb;
          mnx[xnew][0] += abskb;
          mnx[xnew][1] += abskbnew;
        }
      }
    }
  }

}

//_________________________________________________________________________
void nsweeps(int ns)
{
  for(int is=1; is<=ns; is++)
  {
    metropolis_k();
    metropolis_kbar();
    metropolis_plaquette();
#ifndef KAPPA0
    metropolis_s();
    metropolis_sbar();
    metropolis_meson();
#endif
    metropolis_winding_loop();
  }
}

#endif
