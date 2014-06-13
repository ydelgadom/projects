
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator "ranlxd"
!
! See the notes 
!
!    "User's guide for ranlxs and ranlxd [F90 programs]" (December 1997)
!
!    "Double precision implementation of the random number 
!     generator ranlux" (December 1997)
!
! for a detailed description
!
! The publicly accessible subroutines in this module are 
!
! ranlxd(r)
!   real(kind=dp),dimension(:),intent(out) :: r
!   Computation of an array r of double-precision random numbers
!
! rlxd_init(level,seed)
!   integer,intent(in) :: level,seed
!   Initialization of the generator
!
! rlxd_get(state)
!   integer,dimension(:),intent(out) :: state
!   Returns the state of the generator in the form of an integer
!   array with 25 elements 
!
! rlxd_reset(state)
!   integer,dimension(:),intent(in) :: state
!   Resets the state of the generator
!
! Version: 2.1
! Author: Martin Luescher <luscher@mail.desy.de>
! Date: 15.12.1997
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ranlxd_generator

! Stefan Solbrig:
! modifications:
! this implicit none statement, just to be extra careful
implicit none


  public :: ranlxd,rlxd_init,rlxd_get,rlxd_reset

! correct for bug in pgi-compiler
!  public :: scale

  private :: define_constants,error,update

  integer,parameter,private :: dp=selected_real_kind(14)

  integer,save,private :: pr,ir,jr,ir_old,init=0
  integer,dimension(0:11),save,private :: next 

  real(kind=dp),save,private :: zero,one,sbase,base,one_bit,carry
  real(kind=dp),dimension(0:11),save,private :: xdbl

contains
!-----------------------------------------------------------
! correct for bug in pgi-compiler
!  real(kind=dp) function scale(a,b) 
!    real(kind=dp),intent(in) :: a
!    integer,intent(in) :: b
!    scale = a*2.0_dp**b
!  end function scale
!------------------------------------------------------------

  subroutine error(no)

    integer,intent(in) :: no

    select case(no)
      case(0)
        print *, "Error in rlxd_init"
        print *, "Arithmetic on this machine is not suitable for ranlxd"
      case(1)
        print *, "Error in rlxd_init"
        print *, "Bad choice of luxury level (should be 1 or 2)"
      case(2)
        print *, "Error in rlxd_init"
        print *, "Bad choice of seed (should be between 1 and 2^31-1)"
      case(3)
        print *, "Error in rlxd_get"
        print *, "Undefined state or improper argument"
      case(4)
        print *, "Error in rlxd_reset"
        print *, "Arithmetic on this machine is not suitable for ranlxd"
      case(5)
        print *, "Error in rlxd_reset"
        print *, "Unexpected input data"
    end select
    print *, "Program aborted"
    stop

  end subroutine error


  subroutine update()

    integer :: k,kmax
    real(kind=dp) :: y1,y2,y3

    k=0

    do 

      if (ir==0) then
        exit
      end if
      k=k+1

      y1=xdbl(jr)-xdbl(ir)
      y2=y1-carry
      if (y2<zero) then
        carry=one_bit
        y2=y2+one
      else
        carry=zero
      end if
      xdbl(ir)=y2
       
      ir=next(ir)
      jr=next(jr)

    end do

    kmax=pr-12

    do
 
      if (k>kmax) then
        exit
      end if
      k=k+12

      y1=xdbl(7)-xdbl(0)
      y1=y1-carry

      y2=xdbl(8)-xdbl(1)
      if (y1<zero) then
        y2=y2-one_bit
        y1=y1+one
      end if
      xdbl(0)=y1
      y3=xdbl(9)-xdbl(2)
      if (y2<zero) then
        y3=y3-one_bit
        y2=y2+one
      end if
      xdbl(1)=y2
      y1=xdbl(10)-xdbl(3)
      if (y3<zero) then
        y1=y1-one_bit
        y3=y3+one
      end if
      xdbl(2)=y3
      y2=xdbl(11)-xdbl(4)
      if (y1<zero) then
        y2=y2-one_bit
        y1=y1+one
      end if
      xdbl(3)=y1
      y3=xdbl(0)-xdbl(5)
      if (y2<zero) then
        y3=y3-one_bit
        y2=y2+one
      end if
      xdbl(4)=y2
      y1=xdbl(1)-xdbl(6)
      if (y3<zero) then
        y1=y1-one_bit
        y3=y3+one
      end if
      xdbl(5)=y3
      y2=xdbl(2)-xdbl(7)
      if (y1<zero) then
        y2=y2-one_bit
        y1=y1+one
      end if 
      xdbl(6)=y1
      y3=xdbl(3)-xdbl(8)
      if (y2<zero) then
        y3=y3-one_bit
        y2=y2+one
      end if
      xdbl(7)=y2
      y1=xdbl(4)-xdbl(9)
      if (y3<zero) then
        y1=y1-one_bit
        y3=y3+one
      end if
      xdbl(8)=y3
      y2=xdbl(5)-xdbl(10)
      if (y1<zero) then
        y2=y2-one_bit
        y1=y1+one
      end if
      xdbl(9)=y1
      y3=xdbl(6)-xdbl(11)
      if (y2<zero) then
        y3=y3-one_bit
        y2=y2+one
      end if
      xdbl(10)=y2
      if (y3<zero) then
        carry=one_bit
        y3=y3+one
      else
        carry=zero
      end if
      xdbl(11)=y3

    end do

    kmax=pr-k

    do k=1,kmax

      y1=xdbl(jr)-xdbl(ir)
      y2=y1-carry
      if (y2<zero) then
        carry=one_bit
        y2=y2+one
      else
        carry=zero
      end if
      xdbl(ir)=y2
       
      ir=next(ir)
      jr=next(jr)

    end do

    ir_old=ir

  end subroutine update


  subroutine define_constants()

    integer :: k

    init=1
    zero=0.0_dp
    one=1.0_dp

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! this does not work with the pgi-compiler on a
    ! pc running linux. 
    ! pgi- "scale" only works for single precsision
    ! the NAGware f95 compiler works ok
    sbase=scale(one,24)
    base=scale(one,48)
    one_bit=scale(one,-48)

    ! if you have to use pgi use this instead:
    ! sbase=one*2**24
    ! base=one*2**48
    ! one_bit=one*2**(-48)
    ! this is done by redefining the scale function

    do k=0,11

      next(k)=modulo(k+1,12)

    end do

  end subroutine define_constants


  subroutine rlxd_init(level,seed)

    integer,intent(in) :: level,seed
    integer,dimension(0:30) :: xbit
    integer :: ibit,jbit,i,k,l
    real(kind=dp) :: x,y 

    if ((huge(1)<2147483647).or.(radix(1.0_dp)/=2).or. &
    (digits(1.0_dp)<48)) then
      call error(0)
    end if

    select case(level)
      case(1)
        pr=202
      case(2)
        pr=397
      case default
        call error(1)
    end select

    call define_constants()

    i=seed

    do k=0,30

      xbit(k)=modulo(i,2)
      i=i/2
           
    end do

    if ((seed<=0).or.(i/=0)) then
      call error(2)
    end if

    ibit=0
    jbit=18
 
    do k=0,11

      x=zero

      do l=1,48

        y=real(modulo(xbit(ibit)+1,2),dp)
        x=x+x+y

        xbit(ibit)=modulo(xbit(ibit)+xbit(jbit),2)
        ibit=modulo(ibit+1,31)
        jbit=modulo(jbit+1,31)

      end do

      xdbl(k)=one_bit*x

    end do

    carry=zero
    ir=11
    jr=7
    ir_old=0

  end subroutine rlxd_init


  subroutine ranlxd(r)

    real(kind=dp),dimension(:),intent(out) :: r
    integer :: k,kmin,kmax

    if (init==0) then
      call rlxd_init(1,1)
    end if

    kmin=lbound(r,1)
    kmax=ubound(r,1)

    do k=kmin,kmax

      ir=next(ir)

      if (ir==ir_old) then 
             
        call update()

      end if

      r(k)=xdbl(ir)
           
    end do

  end subroutine ranlxd


  subroutine rlxd_get(state)

  integer,dimension(:),intent(out) :: state

    integer :: k
    real(kind=dp) :: x,y1,y2

    if ((init==0).or. &
    ((lbound(state,1)/=1).or.(ubound(state,1)/=25))) then
      call error(3)
    end if

    do k=0,11

      x=sbase*xdbl(k)
      y2=aint(x,dp)
      y1=sbase*(x-y2)
      state(2*k+1)=int(y1)
      state(2*k+2)=int(y2)
    
    end do

    k=12*pr+ir
    k=12*k+jr
    k=12*k+ir_old
    state(25)=2*k+int(carry*base)
  
  end subroutine rlxd_get


  subroutine rlxd_reset(state)

    integer,dimension(:),intent(in) :: state
    integer :: k
    real(kind=dp) :: y1,y2

    if ((huge(1)<2147483647).or.(radix(1.0_dp)/=2).or. &
    (digits(1.0_dp)<48)) then
      call error(4)
    end if

    if ((lbound(state,1)/=1).or.(ubound(state,1)/=25)) then
      call error(5)
    end if

    call define_constants()

    do k=1,24

      if ((state(k)>=int(sbase)).or.(state(k)<0)) then
        call error(5)
      end if

    end do

    k=state(25)
    if (k<0) then
      call error(5)
    end if
    carry=one_bit*real(modulo(k,2),dp)
    k=k/2
    ir_old=modulo(k,12)
    k=k/12
    jr=modulo(k,12)
    k=k/12
    ir=modulo(k,12)
    pr=k/12

    if (((pr/=202).and.(pr/=397)).or.(jr/=modulo((ir_old+7),12))) then
      call error(5)
    end if

    do k=0,11

      y1=real(state(2*k+1),dp)
      y2=real(state(2*k+2),dp)
      xdbl(k)=one_bit*(y1+y2*sbase)
    
    end do

  end subroutine rlxd_reset

end module ranlxd_generator



