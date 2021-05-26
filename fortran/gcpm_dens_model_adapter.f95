! This module implements an adapter for GCPM with a dipole magnetic
! field in solar magnetic (SM) coordinates
module gcpm_dens_model_adapter
  use types
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  implicit none

  ! Types for marshalling.  This is required since the user of this adapter
  ! needs to set additional data for this adapter that is not in the 
  ! interface for funcPlasmaParams() used by the raytracer.
  type :: gcpmStateData
     !	itime	integer	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer :: itime(2)
     !		planetary Kp index
     real(kind=DP) :: akp
     ! Tsyganenko parameters
     real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
     real(kind=DP) :: W1, W2, W3, W4, W5, W6
     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer :: use_igrf

     ! Whether to constrain ray tracing to a constant MLT (1) or not (0)
     integer :: fixed_MLT
     real(kind=DP) :: MLT
     
  end type gcpmStateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: gcpmStateDataP 
     type(gcpmStateData), pointer :: p
  end type gcpmStateDataP

  ! Imported from geopack
  real(kind=SP) :: PSI
  COMMON /GEOPACK1/ PSI

contains

  ! Implementation of the plasma parameters function.
  ! Inputs:
  !   x - position vector in cartesian (SM) coordinates
  ! Outputs:
  !  qs - vector of species charges
  !  Ns - vector of species densities in m^-3
  !  ms - vector of species masses in kg
  ! nus - vector of species collisions in s^-1
  !  B0 - cartesian (SM) background magnetic field in Tesla
  ! In/out:
  ! funcPlasmaParamsData - arbitrary callback data 
  subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
    implicit none

    real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
    real(kind=DP) :: B0(3), B0tmp(3), B0tmp2(3)
    character :: funcPlasmaParamsData(:)

    real(kind=DP) :: p_sm(3)
    real(kind=DP) :: ce,ch,che,co
    real(kind=DP) :: x(3),x_gsm(3)
    
    real(kind=DP) :: r, amlt, alatr
    real(kind=DP) :: outn(4)

    integer :: year, day, hour, min, sec

    real(kind=DP) :: parmod(10)

    integer :: iopt
    real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
    real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE



    type(gcpmStateDataP) :: datap

    iopt = 0

    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

!!$    print *, 'itime=', datap%p%itime
!!$    print *, 'akp=', datap%p%akp
!!$    print *, 'Pdyn=', datap%p%Pdyn
!!$    print *, 'Dst=', datap%p%Dst
!!$    print *, 'ByIMF=', datap%p%BzIMF
!!$    print *, 'W1=', datap%p%W1
!!$    print *, 'W2=', datap%p%W2
!!$    print *, 'W3=', datap%p%W3
!!$    print *, 'W4=', datap%p%W4
!!$    print *, 'W5=', datap%p%W5
!!$    print *, 'W6=', datap%p%W6
!!$    print *, 'use_tsyganenko=', datap%p%use_tsyganenko
!!$    print *, 'use_igrf=', datap%p%use_igrf
!!$    stop
    



    ! Allocate data if not already
    if (.not.(allocated(qs))) then
       allocate(qs(4))
    end if
    if (.not.(allocated(Ns))) then
       allocate(Ns(4))
    end if
    if (.not.(allocated(ms))) then
       allocate(ms(4))
    end if
    if (.not.(allocated(nus))) then
       allocate(nus(4))
    end if
    
    ! Convert from SM x,y,z to GSM x,y,z needed by 
    ! the Tsyganenko model
    call SM_TO_GSM_d(datap%p%itime,x,x_gsm)

    ! convert to spherical coordinates
    p_sm = cartesian_to_spherical(x)
    ! Convert magnetic r,theta,phi to SM local time in hours

    if( datap%p%fixed_MLT == 1 ) then
      amlt = datap%p%MLT
    else
      amlt = mod(24.0_DP*p_sm(2)/(2.0_DP*pi)+12.0_DP,24.0_DP)
    end if
    ! amlt = 0.0_DP
    ! compute r, the geocentric distance in RE
    r = p_sm(1)/R_E
    ! Convert magnetic r,theta,phi to SM latitude in radians
    alatr = pi/2.0_DP - p_sm(3)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !	Input parameters:
    !
    !	itime	integer	dimensions=2
    !		(1) = yearday, e.g. 2001093
    !		(2) = miliseconds of day
    !	r		real(kind=SP)		dimension=1
    !		geocentric radial distance in Re
    !	amlt	real(kind=SP)		dimension=1
    !		solar magnetic local time in hours
    !	alatr	real(kind=SP)		dimension=1
    !		solar magnetic latitude in radians
    !	akp		real(kind=SP)		dimension=1
    !		planetary Kp index
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !	Output parameters:
    !	outn	real(kind=SP)		dimensions=4
    !			(1) = total electron density in 1/cm^3
    !			(2) = total hydrogen density in 1/cm^3
    !			(3) = total helium density in 1/cm^3
    !			(4) = total oxygen density in 1/cm^3
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! print *, 'PRE: itime=', datap%p%itime, 'r=', r, 'amlt=', amlt, 'alatr=', alatr,'akp=',datap%p%akp, 'outn=',outn
    call gcpm_v24(datap%p%itime,r,amlt,&
         alatr,datap%p%akp,outn)
    ! print *, 'POST: itime=', datap%p%itime, 'r=', r, 'amlt=', amlt, 'alatr=', alatr,'akp=',datap%p%akp, 'outn=',outn

    ! print *, 'we did it'
    ce = dble(outn(1))
    ch = dble(outn(2))
    che =dble(outn(3))
    co = dble(outn(4))
    
    if( ce == 0.0_DP ) then
       ce = 1e-12
    end if
    if( ch == 0.0_DP ) then
       ch = 1e-12
    end if
    if( che == 0.0_DP ) then
       che = 1e-12
    end if
    if( co == 0.0_DP ) then
       co = 1e-12
    end if
   
    qs = 1.602e-19_DP*(/ -1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP /);
    ms = (/ 9.10938188e-31_DP, 1.6726e-27_DP, &
         4.0_DP*1.6726e-27_DP, 16.0_DP*1.6726e-27_DP /);
    ! Convert to m^-3;
    Ns = 1.0e6_DP*(/ ce, ch, che, co /);
    nus = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /);

    ! Tsyganenko magnetic field
    
    ! Convert the itime parameter into the Tsyganenko parameters
    year = datap%p%itime(1)/1000
    day = mod(datap%p%itime(1),1000)
    hour = datap%p%itime(2)/(1000*60*60)
    min = (datap%p%itime(2)-hour*(1000*60*60))/(1000*60)
    sec = (datap%p%itime(2)-hour*(1000*60*60)-min*(1000*60))/(1000)

    ! Set the Tsyganenko parameters
    parmod(1) = datap%p%Pdyn   !Pdyn:  between 0.5 and 10 nPa,
    parmod(2) = datap%p%Dst    !Dst:  between -100 and +20,
    parmod(3) = datap%p%ByIMF  !ByIMF: between -10 and +10 nT.
    parmod(4) = datap%p%BzIMF  !BzIMF: between -10 and +10 nT.
    parmod(5) = datap%p%W1     !
    parmod(6) = datap%p%W2     !
    parmod(7) = datap%p%W3     !
    parmod(8) = datap%p%W4     !
    parmod(9) = datap%p%W5     !
    parmod(10)= datap%p%W6     !

    ! Necessary call for the Tsyganenko geopack tools.  Also updates
    ! the common variable psi
    call tsy_recalc(year, day, hour, min, sec)
    if( datap%p%use_igrf == 1 ) then
       ! Find IGRF components in GSM coordinates
       call IGRF_GSM (&
            real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
            B0xBASE,B0yBASE,B0zBASE)
    else
       ! Find the dipole field in SM coordinates
       B0tmp = bmodel_cartesian( x )
       ! Rotate to GSM and convert to nT for convenience with below
       call SM_TO_GSM_d(datap%p%itime,B0tmp,B0tmp2)
       B0xBASE = real(1.0e9_DP*B0tmp2(1))
       B0yBASE = real(1.0e9_DP*B0tmp2(2))
       B0zBASE = real(1.0e9_DP*B0tmp2(3))
    end if
    if( datap%p%use_tsyganenko == 1 ) then
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
    call GSM_TO_SM_d(datap%p%itime,B0tmp,B0)

    
  end subroutine funcPlasmaParams


end module gcpm_dens_model_adapter
