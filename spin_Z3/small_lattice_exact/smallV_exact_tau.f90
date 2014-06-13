
!------------------------------------------------------------------------

   ! 
   !  Program smallV_exact.f90
   !
   
   !  Implements an exact evaluation of the effective center model  
   !  on a 2 x 2 x 3 lattice.
   
   !  C. Gattringer, Nov. 2011 
   
   !  gfortran -O3 smallV_exact.f90  

!------------------------------------------------------------------------
 
  program smallV_exact
  
  implicit none
  
  integer :: s111,s112,s113,s121,s122,s123,  &
             s211,s212,s213,s221,s222,s223
	    
  real :: tau,kappa,mu,eta,etabar
  real :: smu,cmu,snn,smc,sms
  complex :: zz,ss,pp,boltz,ssum,s2sum,psum,p2sum	   	     
  real :: costab(-2:2),sintab(-1:1)
  real :: uu,cc,mm,chi
  integer :: i
  	     
  print*
  print*, "Program smallV_exact.f90 "
  print*
  print*
  print*, "Lattice size = 2 x 2 x 3 "
  print*
  print*, "kappa = "
  read*, kappa
  print*, "mu    = "
  read*, mu
  print*

  costab(-2) = -0.5
  costab(-1) = -0.5
  costab( 0) =  1.0
  costab( 1) = -0.5
  costab( 2) = -0.5
  sintab(-1) = -0.866025403784439
  sintab( 0) =  0.0 
  sintab( 1) =  0.866025403784439

  cmu = kappa*cosh(mu)
  smu = kappa*sinh(mu)
  
  open(unit=14,file="smallV_exact_tau.dat",status="unknown")    
  
  do i = 0,300
    
  zz  = cmplx(0.0,0.0)
  ssum  = cmplx(0.0,0.0)
  s2sum = cmplx(0.0,0.0)
  psum  = cmplx(0.0,0.0)
  p2sum = cmplx(0.0,0.0)

  tau = 0.001*i
  
  do s111 = -1,1
  do s112 = -1,1
  do s113 = -1,1
  do s121 = -1,1
  do s122 = -1,1
  do s123 = -1,1
  do s211 = -1,1
  do s212 = -1,1
  do s213 = -1,1
  do s221 = -1,1
  do s222 = -1,1
  do s223 = -1,1

  snn = costab(s111-s112)+costab(s111-s121)+costab(s111-s211)    &
      + costab(s112-s113)+costab(s112-s122)+costab(s112-s212)    &
      + costab(s113-s111)+costab(s113-s123)+costab(s113-s213)    &
      + costab(s121-s122)+costab(s121-s111)+costab(s121-s221)    &
      + costab(s122-s123)+costab(s122-s112)+costab(s122-s222)    &
      + costab(s123-s121)+costab(s123-s113)+costab(s123-s223)    &
      + costab(s211-s212)+costab(s211-s221)+costab(s211-s111)    &
      + costab(s212-s213)+costab(s212-s222)+costab(s212-s112)    &
      + costab(s213-s211)+costab(s213-s223)+costab(s213-s113)    &
      + costab(s221-s222)+costab(s221-s211)+costab(s221-s121)    &
      + costab(s222-s223)+costab(s222-s212)+costab(s222-s122)    &
      + costab(s223-s221)+costab(s223-s213)+costab(s223-s123)   
	      
  smc = costab(s111)+costab(s112)+costab(s113)+costab(s121)      &
      + costab(s122)+costab(s123)+costab(s211)+costab(s212)      &
      + costab(s213)+costab(s221)+costab(s222)+costab(s223) 
      
  sms = sintab(s111)+sintab(s112)+sintab(s113)+sintab(s121)      &
      + sintab(s122)+sintab(s123)+sintab(s211)+sintab(s212)      &
      + sintab(s213)+sintab(s221)+sintab(s222)+sintab(s223)      
	      
  ss = cmplx(2*tau*snn+2*cmu*smc,2*smu*sms)	
  pp = cmplx(smc,sms)
  
  boltz = exp(ss)
  
  zz    = zz + boltz
  ssum  = ssum - ss*boltz
  s2sum = s2sum + ss**2*boltz
  psum  = psum + pp*boltz !psum + abs(pp)*boltz
  p2sum = p2sum + pp**2*boltz !abs(pp)**2*boltz      
	      
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  
  uu  = - ssum/(12*zz)
  cc  = s2sum/(12*zz) - (ssum/zz)**2/12
  mm  = psum/(12*zz)
  chi = p2sum/(12*zz) - (psum/zz)**2/12 
  
  write(14,100) tau,uu,cc,mm,chi 
  
  enddo

  close(14)
  
  100 format(f8.4,4(2x,f12.6))

  print*
  print*, "Done."
  print*
  
  end program smallV_exact

!------------------------------------------------------------------------
 
