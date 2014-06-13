
!------------------------------------------------------------------------

   ! 
   !  Program local_metro.f90
   !
   
   !  Implements a simulation of the effective center model  
   !  using a worm algorithm. Version for unformatted output.
   
   !  C. Gattringer, July 2010 
   
   !  gfortran -O3 local_metro.f90  

!------------------------------------------------------------------------
 
   !
   !  Modules:
   !

!------------------------------------------------------------------------

   !
   !  Modules for the ranlux random generator.
   !
   
  include "ranlxd_generator.f90"

!------------------------------------------------------------------------

  module lattice

   !
   !  Contains the lattice size and fields for neighbor information.
   !
 
  implicit none
  save

  integer, parameter :: leng = 6  
  integer, parameter :: nsite = leng**3
  integer, dimension(6,nsite) :: neib

  end module lattice
  
!------------------------------------------------------------------------

  module runparams
 
  implicit none
  save

   !
   !  Module for the run parameters.
   !

  integer :: nequi, nmeas, nskip, iseed
  real :: tau,tau0,dtau 
  real :: kappa
  real :: costab(-2:2), sintab(-2:2)
  integer :: ntau
  character(len=128) :: outfile

  end module runparams

!------------------------------------------------------------------------

  module fields

   !
   !  Contains the dimer field.
   !

  use lattice
  implicit none
  save

  integer, dimension(nsite) :: nn
 
  end module fields
  
 !------------------------------------------------------------------------

  module observables

   !
   !  Observables.
   !
 
  implicit none
  save
 
  real, allocatable, dimension(:) :: ee,pp 
 
  end module observables

!------------------------------------------------------------------------

   !
   !  Main program:
   !

!------------------------------------------------------------------------

  program local_metro
  
  use runparams 
  use lattice
  use observables 
  use fields 
  implicit none
  
  integer :: itau,imeas

  print*
  print*, "Program local_metro.f90 "
  print*
  write(*,'(a,i3)') " Compiled with leng =", leng
  print*
  print*

  call readparams  
  call init 
  call initranlux
 
     
  allocate(ee(nmeas)) 
  allocate(pp(nmeas)) 
  
  nn = 0
     
  do itau = 0,ntau  
    
    tau = tau0 + itau*dtau  
     
    write(*,'(a,f8.5)') " tau = ", tau
    
    call donsweeps(nequi)
  
    do imeas = 1,nmeas
      call donsweeps(nskip)
      call measure(imeas)
    enddo  
    
    write(2) ee
    write(2) pp 
  
  enddo
  
  close(2)
  
  deallocate(ee)
  deallocate(pp)
   
  print*
  print*, "Done."
  print*
  
  end program local_metro

!------------------------------------------------------------------------

   !
   !  Subroutines:
   !
   
!------------------------------------------------------------------------

  subroutine measure(im)
  
   !
   !  Computes observables.
   !
  
  use lattice
  use observables
  use runparams
  use fields
  implicit none
  
  integer, intent(in) :: im
  real :: esum,psum
  integer :: is
  
  esum = 0.0
  psum = 0.0
  
  do is = 1,nsite
  
    esum = esum - 2*tau*( costab(nn(is) - nn(neib(1,is)))    &
                        + costab(nn(is) - nn(neib(2,is)))    &
                        + costab(nn(is) - nn(neib(3,is))) )  &
                - 2*kappa*costab(nn(is)) 
    psum = psum + costab(nn(is)) 
    
  enddo
  
  ee(im) = esum
  pp(im) = psum
   
  end subroutine measure
       
!------------------------------------------------------------------------

  subroutine donsweeps(nsweeps)
  
   !
   !  Does nsweeps sweeps.
   !
  
  use lattice
  use ranlxd_generator
  use runparams
  use fields
  implicit none
  
  integer, parameter :: dp=selected_real_kind(14)
  integer, intent(in) :: nsweeps
  integer :: isweeps,is,nnnew
  real :: delta
  real(kind=dp), dimension(nsite) :: ran1,ran2
    
  do isweeps = 1,nsweeps
  
    call ranlxd(ran1)
    call ranlxd(ran2)
    
    do is = 1,nsite
    
      nnnew = floor(ran1(is)*3) - 1
      
      if (nnnew .ne. nn(is) ) then
      
        delta = 2*tau*( costab(nnnew - nn(neib(1,is)))     &
                      + costab(nnnew - nn(neib(2,is)))     &
                      + costab(nnnew - nn(neib(3,is)))     &
		      + costab(nnnew - nn(neib(4,is)))     &
                      + costab(nnnew - nn(neib(5,is)))     &
                      + costab(nnnew - nn(neib(6,is)))     &
		      - costab(nn(is) - nn(neib(1,is)))    &
                      - costab(nn(is) - nn(neib(2,is)))    &
                      - costab(nn(is) - nn(neib(3,is)))    &
		      - costab(nn(is) - nn(neib(4,is)))    &
                      - costab(nn(is) - nn(neib(5,is)))    &
                      - costab(nn(is) - nn(neib(6,is))) )  &
              + 2*kappa*( costab(nnnew) - costab(nn(is)) )
		
	if ( ran2(is) .lt. exp(delta) ) nn(is) = nnnew	
	
      endif		
 
    enddo
     
  enddo
    
  end subroutine donsweeps
    
 !------------------------------------------------------------------------

  subroutine readparams
  
   !
   !  Reads the parameters of the run.
   !
  
  use lattice
  use runparams
  use observables
  implicit none
    
  open(unit = 1, file = "local_metro.start", status = "old")
  
  read(1,*)
  read(1,*) tau0
  read(1,*)
  read(1,*) dtau
  read(1,*)
  read(1,*) ntau
  read(1,*)
  read(1,*) kappa
  read(1,*)
  read(1,*) nequi
  read(1,*)
  read(1,*) nmeas
  read(1,*)
  read(1,*) nskip
  read(1,*)
  read(1,*) iseed
  read(1,*)
  read(1,*) outfile
  
  close(1)
  
  print*
  write(*,'(a,f16.8)')  " tau0    = ", tau0
  write(*,'(a,f16.8)')  " dtau    = ", dtau
  write(*,'(a,i10)')    " ntau    = ", ntau  
  write(*,'(a,f16.8)')  " kappa   = ", kappa
  write(*,'(a,i10)')    " nequi   = ", nequi
  write(*,'(a,i10)')    " nmeas   = ", nmeas
  write(*,'(a,i10)')    " nskip   = ", nskip
  write(*,'(a,i10)')    " iseed   = ", iseed
  write(*,'(a,a)')      " outfile = ", outfile
  print*

   
  open(unit = 2, file = outfile, status = "unknown", form = "unformatted")
 
  write(2)
  write(2) "leng"
  write(2) leng
  write(2) "tau0"
  write(2) tau0
  write(2) "dtau"
  write(2) dtau
  write(2) "ntau"
  write(2) ntau
  write(2) "kappa"
  write(2) kappa
  write(2) "nequi"
  write(2) nequi
  write(2) "nmeas"
  write(2) nmeas
  write(2) "nskip"
  write(2) nskip
  write(2) "iseed"
  write(2) iseed
  write(2)
   
  end subroutine readparams
   
!------------------------------------------------------------------------

  subroutine init

   !
   !  Initializes the neighborhood array and the cosine table.
   !

  use lattice
  use runparams
  implicit none

  integer :: i1,i2,i3,i1p,i2p,i3p,i1m,i2m,i3m,     &
             is,isp1,isp2,isp3,ism1,ism2,ism3


  do i1 = 1,leng
    i1p = i1 + 1
    i1m = i1 - 1
    if (i1 .eq. leng) i1p = 1
    if (i1 .eq. 1) i1m = leng

  do i2 = 1,leng
    i2p = i2 + 1
    i2m = i2 - 1
    if (i2 .eq. leng) i2p = 1
    if (i2 .eq. 1) i2m = leng

  do i3 = 1,leng
    i3p = i3 + 1
    i3m = i3 - 1
    if (i3 .eq. leng) i3p = 1
    if (i3 .eq. 1) i3m = leng
 
  ! compute the site address and the addresses of the sites shifted
  ! by one unit in each direction

  is = i1 + (i2-1)*leng + (i3-1)*leng**2  

  isp1 = i1p + (i2-1)*leng + (i3-1)*leng**2   
  isp2 = i1 + (i2p-1)*leng + (i3-1)*leng**2  
  isp3 = i1 + (i2-1)*leng + (i3p-1)*leng**2  
 
  ism1 = i1m + (i2-1)*leng + (i3-1)*leng**2  
  ism2 = i1 + (i2m-1)*leng + (i3-1)*leng**2  
  ism3 = i1 + (i2-1)*leng + (i3m-1)*leng**2  
 

  ! fill in the neighborhood array

  neib(1,is) = isp1
  neib(2,is) = isp2
  neib(3,is) = isp3
 
  neib(4,is) = ism1
  neib(5,is) = ism2
  neib(6,is) = ism3
 
  end do
  end do
  end do
  
  costab(-2) = -0.5
  costab(-1) = -0.5
  costab( 0) =  1.0
  costab( 1) = -0.5
  costab( 2) = -0.5

  costab(-2) = -0.5
  costab(-1) = -0.5
  costab( 0) =  1.0
  costab( 1) = -0.5
  costab( 2) = -0.5
    
  print*, "Neighbor array and cosine table initialized."

  end subroutine init

!------------------------------------------------------------------------

  subroutine initranlux

   !
   !  Tests parameters of the ranlux random generator and initializes it.
   !

  use runparams
  use ranlxd_generator
  implicit none
  
  integer, parameter :: dp=selected_real_kind(14)

  real(kind=dp), parameter :: theOne = 1.0_dp
  real(kind=dp) :: theBig, theSmall
 
  theBig=scale(theOne,48)
  theSmall=scale(theOne,-48)
  
  ! print *,theOne,theBig,theSmall                                      
  
  if( theBig == 0.0_dp .or. theBig == 1.0_dp .or. &
      theSmall == 0.0_dp .or. theSmall == 1.0_dp ) then
    print *," intrinsic funtion SCALE is not working correctly "
    print *," include ranlxd_generator-pgi.f90 instead         "
    STOP
   
  end if
    
  call rlxd_init(2,iseed)
  
  print*, "Parameters of ranlux tested and ranlux initialized."
  print*
  
  end subroutine initranlux

!------------------------------------------------------------------------
