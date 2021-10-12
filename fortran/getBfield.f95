! This program dumps bfield at locations
program getBfield
  use types
  use util
  use constants, only : R_E, PI, VERSION
  use bmodel_dipole
  USE ISO_FORTRAN_ENV ! for OUTPUT_UNIT definition, from 2003 standard
  implicit none

  ! set up for in/out parameters
  character(len=10000) :: infilename
  character(len=10000) :: filename
  integer,parameter :: outfile=50
  integer,parameter :: infile=50
  character (len=10000) :: buffer
  real(kind=DP) :: tmpinput
  integer :: foundopt, k, nrays, status

  ! Whether to use (1) or not use (0) the tsyganenko corrections
  integer :: use_tsyganenko
  ! Whether to use (1) IGRF or not use (0) and use dipole instead
  integer :: use_igrf

  !	itime	integer	dimensions=2
  !		(1) = yearday, e.g. 2001093
  !		(2) = miliseconds of day
  integer :: itime(2)
  integer :: year, day, hour, min, sec
  ! Tsyganenko parameters
  real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
  real(kind=DP) :: W1, W2, W3, W4, W5, W6
  real(kind=DP) :: parmod(10)
  real(kind=SP) :: psi
  integer :: iopt
  ! Tsyganenko corrections
  real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
  ! Base B field 
  real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE
  ! input positions
  real(kind=DP) :: x(3),x_gsm(3),pos(3)
  ! B output
  real(kind=DP) :: B0(3), B0tmp(3), B0tmp2(3)

  ! final array -- f
  real(kind=DP), allocatable :: f(:,:)
  real(kind=DP), allocatable :: pos_a(:,:)

  ! ---------------------------------- get inputs from input arguments ----------------------------------
  call getopt_named( 'filename', filename, foundopt )
  call getopt_named( 'infilename', infilename, foundopt )

  ! yearday
  call getopt_named( 'yearday', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     itime(1) = tmpinput
  end if
  ! milliseconds_day
  call getopt_named( 'milliseconds_day', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     itime(2) = tmpinput
  end if
  ! use_tsyganenko
  call getopt_named( 'use_tsyganenko', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     use_tsyganenko = floor(tmpinput)
  end if
  ! use_igrf
  call getopt_named( 'use_igrf', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     use_igrf = floor(tmpinput)
  end if
  ! tsyganenko_Pdyn
  call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     Pdyn = tmpinput
  end if
  ! tsyganenko_Dst
  call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     Dst = tmpinput
  end if
  ! tsyganenko_ByIMF
  call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     ByIMF = tmpinput
  end if
  ! tsyganenko_BzIMF
  call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     BzIMF = tmpinput
  end if
  ! tsyganenko_W1
  call getopt_named( 'tsyganenko_W1', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W1 = tmpinput
  end if
  ! tsyganenko_W2
  call getopt_named( 'tsyganenko_W2', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W2 = tmpinput
  end if
  ! tsyganenko_W3
  call getopt_named( 'tsyganenko_W3', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W3 = tmpinput
  end if  
  ! tsyganenko_W4
  call getopt_named( 'tsyganenko_W4', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W4 = tmpinput
  end if
  ! tsyganenko_W5
  call getopt_named( 'tsyganenko_W5', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W5 = tmpinput
  end if
  ! tsyganenko_W6
  call getopt_named( 'tsyganenko_W6', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     W6 = tmpinput
  end if
  ! nrays
  call getopt_named( 'nrays', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     nrays = floor(tmpinput)
  end if

  ! allocate input and output - same size
  allocate(pos_a(3, nrays))
  allocate(f(3, nrays))

  ! get input positions
  open(unit=infile, file=infilename, status="old")

  do k=1,nrays
     ! Read in ray positions
     read(infile, *,iostat=status), pos
     if( status /= 0 ) then
        ! end of file or error, abort the loop
        exit
     end if
     ! allocate to position
     pos_a(:,k) = pos
  end do
  flush(OUTPUT_UNIT)
     
  ! ---------------------------------- get Bfield ----------------------------------
  do k=1,nrays
     flush(OUTPUT_UNIT)

     ! get input position
     x(1) = pos_a(1,k)
     x(2) = pos_a(2,k)
     x(3) = pos_a(3,k)

     ! Convert from SM x,y,z to GSM x,y,z needed by 
     ! the Tsyganenko model
     call SM_TO_GSM_d(itime,x,x_gsm)

     ! Convert the itime parameter into the Tsyganenko parameters
     year = itime(1)/1000
     day = mod(itime(1),1000)
     hour = itime(2)/(1000*60*60)
     min = (itime(2)-hour*(1000*60*60))/(1000*60)
     sec = (itime(2)-hour*(1000*60*60)-min*(1000*60))/(1000)

     ! Set the Tsyganenko parameters
     parmod(1) = Pdyn   !Pdyn:  between 0.5 and 10 nPa,
     parmod(2) = Dst    !Dst:  between -100 and +20,
     parmod(3) = ByIMF  !ByIMF: between -10 and +10 nT.
     parmod(4) = BzIMF  !BzIMF: between -10 and +10 nT.
     parmod(5) = W1     !
     parmod(6) = W2     !
     parmod(7) = W3     !
     parmod(8) = W4     !
     parmod(9) = W5     !
     parmod(10)= W6     !
   
     ! necessary call
     call tsy_recalc(year, day, hour, min, sec)

     if( use_igrf == 1 ) then
        ! Find IGRF components in GSM coordinates
        call IGRF_GSM (&
         real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
         B0xBASE,B0yBASE,B0zBASE)
     else
        ! Find the dipole field in SM coordinates
        B0tmp = bmodel_cartesian( x )
        ! Rotate to GSM and convert to nT for convenience with below
        call SM_TO_GSM_d(itime,B0tmp,B0tmp2)
        B0xBASE = real(1.0e9_DP*B0tmp2(1))
        B0yBASE = real(1.0e9_DP*B0tmp2(2))
        B0zBASE = real(1.0e9_DP*B0tmp2(3))
     end if

     ! Tsyganenko magnetic field
     if( use_tsyganenko == 1 ) then
        call T04_s( iopt, real(parmod), real(psi), &
          real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
          B0xTsy, B0yTsy, B0zTsy)
     else
        B0xTsy = 0.0
        B0yTsy = 0.0
        B0zTsy = 0.0
     end if
    
   ! Add the field and Tsyganenko corrections together and convert from
   ! nT to T
   B0tmp(1) = (B0xBASE+B0xTsy)*1.0e-9_DP
   B0tmp(2) = (B0yBASE+B0yTsy)*1.0e-9_DP
   B0tmp(3) = (B0zBASE+B0zTsy)*1.0e-9_DP

   ! We're in GSM coordinates.  Rotate back to SM
   call GSM_TO_SM_d(itime,B0tmp,B0)

   ! save output to f
   f(:,k) = B0
   ! print *, B0
  end do

  ! ---------------------------------- write to out file ----------------------------------
  open(unit=outfile, file=filename, status="replace")
  ! Print out the data
  write(outfile, '(es24.15e3)'), f
  
  flush(infile)
  flush(outfile)
  ! close
  close(unit=infile)
  close(unit=outfile)
  
end program getBfield