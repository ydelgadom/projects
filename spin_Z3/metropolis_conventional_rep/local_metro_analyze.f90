
!------------------------------------------------------------------------

   ! 
   !  Program local_metro_analyze.f90
   !
   
   !  Anlyzes results from the metropolis simulation of the effective center 
   !  model in the full case. Version for unformatted output.  
   
   !  C. Gattringer, July 2010 
   
   !  gfortran -O3 local_metro_analyze.f90  


!------------------------------------------------------------------------
 
   !
   !  Modules:
   !

!------------------------------------------------------------------------
   
  module obsdata

   !
   !  Contains the raw data and parameters.
   !
 
  implicit none
  save

  integer :: nsite,ntau,nmeas
  real, allocatable, dimension(:) :: tau
  real, allocatable, dimension(:,:) :: ee,pp 
 
  end module obsdata
  
!------------------------------------------------------------------------
  
  module finaldata

   !
   !  Contains the raw data and parameters.
   !
 
  use obsdata
  implicit none
  save

  integer, parameter :: ntaumax = 12

  real, dimension(0:ntaumax) :: eaver,eerr
  real, dimension(0:ntaumax) :: caver,cerr
  real, dimension(0:ntaumax) :: paver,perr  
  real, dimension(0:ntaumax) :: saver,serr
    
  end module finaldata
  
!------------------------------------------------------------------------

   !
   !  Main program:
   !

!------------------------------------------------------------------------

  program local_metro_analyze

  use obsdata
  use finaldata   
  implicit none
  
  integer :: itau

  print*
  print*, "Program local_metro_analyze.f90 "
  print*

  call readdata
  call finalize
  
  open(unit = 1, file = "local_metro_analyze.outobs", status = "unknown")

  do itau = 0,ntau
  
  write(1,*) tau(itau), eaver(itau), eerr(itau),  &
                        caver(itau), cerr(itau),  &
                        paver(itau), perr(itau),  &  
                        saver(itau), serr(itau) 
  enddo			  
 
  write(1,*)
  
  close(1)
  
  100 format('(f8.5,4(4x,e16.6,1x,e16.6))')		  

  print*
  print*, "Done."
  print*

  end program local_metro_analyze
  
!------------------------------------------------------------------------

   !
   !  Subroutines:
   !

!------------------------------------------------------------------------

  subroutine finalize
  
   !
   !  Finalizes the data.
   !
  
  use obsdata
  use finaldata
  implicit none
  
  integer :: it,im
  real :: esum,psum,efsum,pfsum,efblock,pfblock,cfsum,sfsum 
  
  do it = 0,ntau
  
    esum = 0.0
    psum = 0.0
 
    do im = 1,nmeas

      esum = esum + ee(it,im)
      psum = psum + pp(it,im)
       
    enddo
  
    eaver(it) = esum/nmeas
    paver(it) = psum/nmeas
         
    efsum = 0.0
    pfsum = 0.0
 
    do im = 1,nmeas
  
      efsum = efsum + ( eaver(it) - ee(it,im) )**2
      pfsum = pfsum + ( paver(it) - pp(it,im) )**2 
       
    enddo  
  
    eerr(it) = sqrt(efsum)/nmeas
    perr(it) = sqrt(pfsum)/nmeas
     
    caver(it) = efsum/nmeas  
    saver(it) = pfsum/nmeas
   
    cfsum = 0.0
    sfsum = 0.0
 
    do im = 1,nmeas
  
      efblock = ( efsum - ( eaver(it) - ee(it,im) )**2 )/(nmeas-1)
      pfblock = ( pfsum - ( paver(it) - pp(it,im) )**2 )/(nmeas-1)
       
      cfsum = cfsum + ( caver(it) - efblock )**2
      sfsum = sfsum + ( saver(it) - pfblock )**2
             
    enddo  
  
    cerr(it) = sqrt(cfsum)
    serr(it) = sqrt(sfsum)
   
  enddo  
  
  eaver = eaver/nsite
  eerr  = eerr/nsite
  paver = paver/nsite
  perr  = perr/nsite
  caver = caver/nsite
  cerr  = cerr/nsite  
  saver = saver/nsite
  serr  = serr/nsite
    
  end subroutine finalize
   
!------------------------------------------------------------------------

  subroutine readdata
  
   !
   !  Reads the data.
   !
   
  use obsdata
  implicit none
  
  character(len=128) :: infile
  integer :: leng,nequi,nskip,iseed,itau,imeas
  real :: tau0,dtau,kappa
  
  print*
  print*, "Input file ="
  read*, infile
  print*
  
  open(unit = 1, file = infile, status = "old", form = "unformatted")
  
  read(1)  
  read(1)  
  read(1) leng 
  read(1)  
  read(1) tau0
  read(1)
  read(1) dtau
  read(1)
  read(1) ntau
  read(1)   
  read(1) kappa
  read(1)
  read(1) nequi
  read(1)
  read(1) nmeas
  read(1)
  read(1) nskip
  read(1)
  read(1) iseed
  read(1)
   
  write(*,'(a,i10)')    " leng  = ", leng
  write(*,'(a,f16.8)')  " tau0  = ", tau0
  write(*,'(a,f16.8)')  " dtau  = ", dtau
  write(*,'(a,i10)')    " ntau  = ", ntau 
  write(*,'(a,f16.8)')  " kappa = ", kappa
  write(*,'(a,i10)')    " nequi = ", nequi
  write(*,'(a,i10)')    " nmeas = ", nmeas
  write(*,'(a,i10)')    " nskip = ", nskip
  write(*,'(a,i10)')    " iseed = ", iseed
  print*
  
  nsite = leng**3
     
  allocate(ee(0:ntau,nmeas))  
  allocate(pp(0:ntau,nmeas))   		
  allocate(tau(0:ntau))   		
  		   
  do itau = 0,ntau  
    
    tau(itau) = tau0 + itau*dtau  

    print*, "reading tau = ", tau(itau)

    read(1) ee(itau,:)
    read(1) pp(itau,:)
     
  enddo
  
  close(1)
    
  end subroutine readdata 
  
!------------------------------------------------------------------------

 
