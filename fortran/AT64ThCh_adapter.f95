! This is a density adapter for the AT64ThCh model
module AT64ThCh_adapter
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  implicit none

  ! Types for marshalling.  This is required since the user of this
  ! adapter often needs to set additional state or configuration data
  ! for that is not in the interface for funcPlasmaParams() used by
  ! the raytracer.
  !
  ! You can put whatever data you like in here.  Included below are
  ! some parameters for the Tsyganenko model.  These should be set by
  ! the caller, associated with a pointer container, also defined
  ! below, then the pointer container type copied into the data(:)
  ! parameter with the TRANSFER intrinsic prior to starting the
  ! raytracer.
  !
  ! For example, in the caller, create an instance of the StateData
  ! type and the pointer container (which will actually be
  ! marshalled):
  !
  !   type(StateData),target :: my_state_data
  !   type(StateDataP) :: my_state_dataP
  ! 
  ! Also create an allocatable character string which will hold the
  ! marshalled data:
  !
  !   character, allocatable :: data(:)
  ! 
  ! Then set the parameters in the caller
  ! (e.g. my_state_data%use_igrf=1).
  ! 
  ! Finally, associate the pointer and marshall as follows:
  ! 
  !   my_state_dataP%p => my_state_data
  !   sz = size(transfer(my_state_datap, data))
  !   allocate(data(sz))
  !   data = transfer(my_state_dataP, data)
  ! 
  type :: AT64ThChStateData
     !	itime	integer	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer :: itime(2)
     integer :: gcpm_kp
     ! Tsyganenko parameters
     real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
     real(kind=DP) :: W1, W2, W3, W4, W5, W6
     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer :: use_igrf
  end type AT64ThChStateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: AT64ThChStateDataP 
     type(AT64ThChStateData), pointer :: p
  end type AT64ThChStateDataP

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

    real(kind=DP) :: ce,ch,che,co, g0, gh
    real(kind=DP) :: x(3),x_gsm(3)
    real(kind=DP) :: tt, T, mpg, radial_dist, Rp, SN, temp_gradient, tran, z, zg
    real(kind=DP) :: R, L, a, c_p, cos_lat, etrans_dens, h, H0, H1, H3, Lpp, Lw
    real(kind=DP) :: n10, n30, ne, ne_tmp, neutral_temp, nh, noh
    real(kind=DP) :: OH_transition_temp, OH_transition_height, peak_height, r0, R13
    real(kind=DP) :: lat_angle, bmag_base, zbrat, b_oh, bmag
    
    ! Tsyganenko parameters
    integer :: year, day, hour, min, sec
    real(kind=DP) :: parmod(10)
    real(kind=SP) :: psi
    integer :: iopt
    ! Tsyganenko corrections
    real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
    ! Base B field 
    real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE
    real(kind=SP) :: B0xoh, B0yoh, B0zoh
    real(kind=SP) :: XX(500), YY(500), ZZ(500)
    real(kind=SP) :: dsmax, rlim, lmax, m, dir, r0s, err,  XF, YF, ZF

    type(AT64ThChStateDataP) :: datap

    iopt = 0

    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

    ! Allocate data if not already
    if (.not.(allocated(qs))) then
       allocate(qs(3))
    end if
    if (.not.(allocated(Ns))) then
       allocate(Ns(3))
    end if
    if (.not.(allocated(ms))) then
       allocate(ms(3))
    end if
    if (.not.(allocated(nus))) then
       allocate(nus(3))
    end if
    
    ! Convert from SM x,y,z to GSM x,y,z needed by 
    ! the Tsyganenko model
    call SM_TO_GSM_d(datap%p%itime,x,x_gsm)
    
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

    !!!!!!!!!!!  FILL IN YOUR DENSITY MODEL HERE -- 
    ! Fill in qs (charge), ms (masses), Ns (number density in m^-3),
    ! and nus (collision frequency, currently unused), for all species
    ! of interest

    ! think we need radial distance first -- should be in meter
    radial_dist = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))

    ! hope this is in m - height from earth
    h = radial_dist - R_E 
    
    ! 400 km up is the transition height
    OH_transition_height = 400.0e3_DP
    ! okay this is in m
    r0 = R_E + OH_transition_height
    
    ! geocentric ratio
    R = radial_dist / r0
   
    ! we need B now (mag) 
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

    bmag = sqrt(B0xBASE*B0xBASE + B0yBASE*B0yBASE + B0zBASE*B0zBASE)

    ! traccing params
    LMAX = 500
    DIR = 1
    DSMAX = 1
    ERR = 0.0001
    RLIM = 60
    R0s = (OH_transition_height+R_E)/R_E

    ! check to make sure the ray is above the transition height
    if (h > 400e3_DP) then

    ! use geopack tools to trace to northern hemisphere 
    ! output coordinates at fieldline footprint at OH transition height
    CALL TRACE_08 (real(x_gsm(1)/R_E),real(x_gsm(2)/R_E),real(x_gsm(3)/R_E), &
      DIR,DSMAX,ERR,RLIM,real(R0s),IOPT,PARMOD,T04_s,IGRF_GSM,XF,YF,ZF,& 
      XX,YY,ZZ,M,LMAX)
   
    ! find B here
    call IGRF_GSM (real(XF), real(YF), real(ZF), B0xoh, B0yoh, B0zoh)
    b_oh = sqrt(B0xoh*B0xoh + B0yoh*B0yoh + B0zoh*B0zoh)

    zbrat = real(bmag/b_oh)
    ! if it's not, just set the zbrat to 1
    else
    zbrat = 1
    end if

    ! very simple dipole implementation to get Lshell
    lat_angle = atan(x(3),x(1))
    ! maybe in radians but honestly no clue
    cos_lat = cos(lat_angle)
    L = (radial_dist / R_E) / (cos_lat*cos_lat)

    ! okay now we need a
    ! there are in K
    temp_gradient = 800.0_DP
    OH_transition_temp = 750.0_DP
    ! use r0 in Mm
    a = temp_gradient * (r0/1.0e6_DP) / OH_transition_temp-1.0_DP

    ! this is in Mm also ??
    tt = (R* ( 1.0_DP + a ) - a) / R

    ! this is now a ratio - a is in Mm and convert r0 to Mm
    zg = (r0/1.0e6_DP) / a * log(tt)
    peak_height = 300.0e3_DP
    ! this is a ratio, just use meters
    Rp = ( R_E + peak_height ) / r0

    ! a is in Mm so this is also in Mm
    c_p = 1.0_DP / ((Rp * (1.0_DP + a) - a ) * Rp)
    neutral_temp = 1000.0_DP 
    ! mass of a proton times grav acceleration [ kg m / s^2 ] -- varying w alt
    gh = 9.80665_DP
    ! gh = g0*(R_E/radial_dist)
    mpg = 1.6726219e-27_DP * gh
    
    ! neutral 0 scale height in Mm
    H0 = 1.380658e-23_DP * neutral_temp / (16.0_DP * mpg) / 1.0e6_DP
    ! this is in Mm, zg is a ratio but cp and H0 are in Mm
    z = zg + ( c_p * H0 * exp( ( (peak_height - h) / 1.0e6_DP ) / H0) )
    T = OH_transition_temp * tt
    ! all in Mm below
    H1 = 1.380658e-23_DP * OH_transition_temp / mpg / 1.0e6_DP
    H3 = 1.380658e-23_DP * OH_transition_temp / (16.0_DP * mpg) / 1.0e6_DP
    etrans_dens = 2.0e11_DP
    n10 = 0.5_DP * etrans_dens
    n30 = 0.5_DP * etrans_dens

    ne_tmp = sqrt( (etrans_dens * OH_transition_temp ) * zbrat * ((n10 * OH_transition_temp) &
    * exp( -1.0_DP * z / H1) + (n30 * OH_transition_temp) * exp( -1.0_DP * z / H3)) ) / T

    R13 = (n10 / n30) * exp( z * (( H1 - H3) / (H1 * H3) ))
    SN = 124.0_DP * (3.0_DP / L)**(4.0_DP) * 1.0e6_DP
    Lpp = 5.6_DP - (0.46_DP * datap%p%gcpm_kp)
    Lw = 0.14_DP
    tran = 0.5_DP * tanh( 3.4534_DP * (L - Lpp) / Lw ) + 0.5_DP
    ne = (1.0_DP - tran) * ne_tmp + tran * SN
    noh = ne / (1.0_DP + R13)
    nh = ne / (1.0_DP + (1.0_DP / R13) )

    ! electron, o+, hydrogen
    qs = 1.602e-19_DP*(/ -1.0_DP, 1.0_DP, 1.0_DP /);
    ms = (/ 9.10938188e-31_DP, 16.0_DP*1.6726e-27_DP, 1.6726e-27_DP /);
    Ns = (/ ne, noh, nh /);
    nus = (/ 0.0_DP, 0.0_DP, 0.0_DP /);

    !!!!!!!!!!!  END FILL IN YOUR DENSITY MODEL HERE

    ! Tsyganenko magnetic field
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


end module AT64ThCh_adapter
