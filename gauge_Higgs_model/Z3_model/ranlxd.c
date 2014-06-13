
/*******************************************************************************
*
* Random number generator "ranlxd"
*
* See the notes 
*
*   "User's guide for ranlxs and ranlxd [C programs]" (December 1997)
*
*   "Double precision implementation of the random number 
*    generator ranlux" (December 1997)
*
* for a detailed description
*
* The externally accessible functions are 
*
*   void ranlxd(double r[],int n)
*     Computes the next n double-precision random numbers and 
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
* 
*   void rlxd_init(int level,int seed)
*     Initialization of the generator
*
*   void rlxd_get(int state[])
*     Extracts the current state of the generator and stores the 
*     information in the array state[25]
*
*   void rlxd_reset(int state[])
*     Resets the generator to the state defined by the array state[25]
*
* Version: 2.1
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 15.12.1997
*
*******************************************************************************/

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

  static int pr,ir,jr,ir_old,init=0;
  static int next[12]; 

  static double zero,one,sbase,base,one_bit,carry;
  static double xdbl[12];

  static void error(no)
  int no;
  {
    switch(no)
    {
      case 0:
        printf("Error in rlxd_init\n");
        printf("Arithmetic on this machine is not suitable for ranlxd\n");
        break;
      case 1:
        printf("Error in rlxd_init\n");
        printf("Bad choice of luxury level (should be 1 or 2)\n");
        break;
      case 2:
        printf("Error in rlxd_init\n");
        printf("Bad choice of seed (should be between 1 and 2^31-1)\n");
        break;
      case 3:
        printf("Error in rlxd_get\n");
        printf("Undefined state\n");
        break;
      case 4:
        printf("Error in rlxd_reset\n");
        printf("Arithmetic on this machine is not suitable for ranlxd\n");
        break;
      case 5:
        printf("Error in rlxd_reset\n");
        printf("Unexpected input data\n");
        break;
    }         
    printf("Program aborted\n");
    exit(0);
  }
  

  static void update()
  {
    int k,kmax;
    double y1,y2,y3;

    for (k=0;ir>0;++k) 
    {
      y1=xdbl[jr]-xdbl[ir];
      y2=y1-carry;
      if (y2<zero) 
      { 
        carry=one_bit;
        y2+=one;
      }
      else 
        carry=zero;
      xdbl[ir]=y2;
      ir=next[ir];
      jr=next[jr];
    }

    kmax=pr-12;

    for (;k<=kmax;k+=12)
    {
      y1=xdbl[7]-xdbl[0];
      y1-=carry;

      y2=xdbl[8]-xdbl[1];
      if (y1<zero)
      {
        y2-=one_bit;
        y1+=one;
      }
      xdbl[0]=y1;
      y3=xdbl[9]-xdbl[2];
      if (y2<zero)
      {
        y3-=one_bit;
        y2+=one;
      }
      xdbl[1]=y2;
      y1=xdbl[10]-xdbl[3];
      if (y3<zero)
      {
        y1-=one_bit;
        y3+=one;
      }
      xdbl[2]=y3;
      y2=xdbl[11]-xdbl[4];
      if (y1<zero)
      {
        y2-=one_bit;
        y1+=one;
      }
      xdbl[3]=y1;
      y3=xdbl[0]-xdbl[5];
      if (y2<zero)
      {
        y3-=one_bit;
        y2+=one;
      }
      xdbl[4]=y2;
      y1=xdbl[1]-xdbl[6];
      if (y3<zero)
      {
        y1-=one_bit;
        y3+=one;
      }
      xdbl[5]=y3;
      y2=xdbl[2]-xdbl[7];
      if (y1<zero)
      {
        y2-=one_bit;
        y1+=one;
      }
      xdbl[6]=y1;
      y3=xdbl[3]-xdbl[8];
      if (y2<zero)
      {
        y3-=one_bit;
        y2+=one;
      }
      xdbl[7]=y2;
      y1=xdbl[4]-xdbl[9];
      if (y3<zero)
      {
        y1-=one_bit;
        y3+=one;
      }
      xdbl[8]=y3;
      y2=xdbl[5]-xdbl[10];
      if (y1<zero)
      {
        y2-=one_bit;
        y1+=one;
      }
      xdbl[9]=y1;
      y3=xdbl[6]-xdbl[11];
      if (y2<zero)
      {
        y3-=one_bit;
        y2+=one;
      }
      xdbl[10]=y2;
      if (y3<zero)
      {
        carry=one_bit;
        y3+=one;
      }
      else
        carry=zero;
      xdbl[11]=y3;
    }  

    kmax=pr;

    for (;k<kmax;++k) 
    {
      y1=xdbl[jr]-xdbl[ir];
      y2=y1-carry;
      if (y2<zero) 
      { 
        carry=one_bit;
        y2+=one;
      }
      else 
        carry=zero;
      xdbl[ir]=y2;
      ir=next[ir];
      jr=next[jr];
    }
    ir_old=ir;
  }


  static void define_constants()
  {
    int k;

    init=1;
    zero=0.0;
    one=1.0;
    sbase=ldexp(one,24);
    base=ldexp(one,48);
    one_bit=ldexp(one,-48);

    for (k=0;k<12;++k) 
    {
      next[k]=(k+1)%12;
    }
  }


  void rlxd_init(level,seed)
  int level,seed;
  {
    int ibit,jbit,i,k,l,xbit[31];
    double x,y;

    if ((INT_MAX<2147483647)||
    (FLT_RADIX!=2)||(DBL_MANT_DIG<48))
      error(0);

    if      (level==1)
      pr=202;
    else if (level==2)
      pr=397;
    else
      error(1);

    define_constants();
    i=seed;

    for (k=0;k<31;++k) 
    {
      xbit[k]=i%2;
      i/=2;
    }

    if ((seed<=0)||(i!=0))
      error(2);

    ibit=0;
    jbit=18;
 
    for (k=0;k<12;++k) 
    {
      x=zero;

      for (l=1;l<=48;++l) 
      {
        y=(double)((xbit[ibit]+1)%2);
        x+=x+y;
        xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
        ibit=(ibit+1)%31;
        jbit=(jbit+1)%31;
      }
      xdbl[k]=one_bit*x;
    }

    carry=zero;
    ir=11;
    jr=7;
    ir_old=0;
  }


  void ranlxd(r,n)
  double r[];
  int n;
  {
    int k;

    if (init==0)
      rlxd_init(1,1);

    for (k=0;k<n;++k) 
    {
      ir=next[ir];
      if (ir==ir_old)
        update();
      r[k]=xdbl[ir];
    }
  }


  void rlxd_get(state)
  int state[];
  {
    int k;
    double x,y1,y2;

    if (init==0)
      error(3);

    for (k=0;k<12;++k) 
    {
      x=sbase*xdbl[k];
      y1=sbase*modf(x,&y2);
      state[2*k]=(int)y1;
      state[2*k+1]=(int)y2;
    }

    k=12*pr+ir;
    k=12*k+jr;
    k=12*k+ir_old;
    state[24]=2*k+(int)(carry*base);
  }


  void rlxd_reset(state)
  int state[];
  {
    int k;
    double y1,y2;

    if ((INT_MAX<2147483647)||
    (FLT_RADIX!=2)||(DBL_MANT_DIG<48))
      error(4);

    define_constants();

    for (k=0;k<24;++k) 
    {
      if ((state[k]>=(int)sbase)||(state[k]<0))
        error(5);
    }

    k=state[24];
    if (k<0)
      error(5);
    carry=one_bit*(double)(k%2);
    k/=2;
    ir_old=k%12;
    k/=12;
    jr=k%12;
    k/=12;
    ir=k%12;
    pr=k/12;

    if (((pr!=202)&&(pr!=397))||(jr!=((ir_old+7)%12)))
      error(5);

    for (k=0;k<12;++k) 
    {
      y1=(double)state[2*k];
      y2=(double)state[2*k+1];
      xdbl[k]=one_bit*(y1+y2*sbase);
    }
  }




