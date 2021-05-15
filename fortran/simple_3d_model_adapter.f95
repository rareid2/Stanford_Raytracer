! This is an empty implementation of a density adapter.  Fill in the pieces
! to adapt this to your own application.
module simple_3d_model_adapter
  use util
  use constants, only : R_E, PI, D2R, R2D, REkm
  use bmodel_dipole
  use pp_profile_d  ! Plasmapause profile (double precision)
  use switch_d, only : switch ! Hyperbolic tangent transition function
  USE ISO_FORTRAN_ENV ! for OUTPUT_UNIT definition, from 2003 standard
  ! use iri_curves
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
  type :: simpleStateData
     !	itime	integer	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer :: itime(2)
     ! Tsyganenko parameters
     real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
     real(kind=DP) :: W1, W2, W3, W4, W5, W6

     ! Kp
     real(kind=DP) :: kp


     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer :: use_igrf

     ! Whether to constrain ray tracing to a constant MLT (1) or not (0)
     integer :: fixed_MLT
     real(kind=DP) :: MLT

     ! ! Which sections to include:
     ! integer :: do_cap
     ! integer :: do_trough
     ! integer :: do_ps
     ! integer :: do_iono

  end type simpleStateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: simpleStateDataP 
     type(simpleStateData), pointer :: p
  end type simpleStateDataP

  ! Imported from geopack
  real(kind=SP) :: PSI
  COMMON /GEOPACK1/ PSI

  real(kind=DP), parameter :: iono_peak_altitude = 350.0 ! altitude of f2 peak
  real(kind=DP), parameter :: iono_merge_radius  = 10000.0 ! Radius in Km of iono - ps transition
  real(kind=DP), parameter :: altrans = 10.0 ! L-shell transition width into cap region
  real(kind=DP), parameter :: rz12= 0.0  ! 13-month sunspot number (~0 for 2010)
  real(kind=DP), parameter :: f107 = 70.0 ! F10.7 (I'm lazy and haven't looked it up yet)
  integer(kind=DP), parameter :: do_cap = 0
  integer(kind=DP), parameter :: do_trough = 1
  integer(kind=DP), parameter :: do_ps = 1
  integer(kind=DP), parameter :: do_iono = 1
  ! Module globals to avoid recalculating the auroral oval. Sure, why not.
  real(kind=DP) :: oldmlt = -1.0
  real(kind=DP) :: oldkp  = -1.0


contains


  function ne_ps(L, a8, a9, iyear, doy, akp)
    ! The Carpenter + Anderson equatorial plasmasphere density
    !  profile (Gallagher 2000, equation 5).
    real(kind=DP) :: L, a8, a9, akp, doy, iyear
    real(kind=DP) :: h,  x234, rz12, a6, a7
    real(kind=DP) ::  ne_ps, doy_factor

    a6 = -0.79
    a7 = 5.208
    ! rz12= 0.0  ! 13-month sunspot number (~0 for 2010)

    ! ---------- Carpenter / Anderson plasmasphere -------------
    h = (1.0_DP + (L/a8)**(2.0_DP*(a9 - 1)))**(-a9/(a9 - 1.0_DP))
    doy_factor=pi*(doy+9.0_DP)/365.0_DP
    x234=( 0.15_DP*(cos(2.0_DP*doy_factor) - 0.5_DP*cos(4.0_DP*doy_factor))+ &
         &(0.00127_DP*rz12-0.0635_DP) ) * exp(-(L-2.0_DP)/1.5_DP)

    ne_ps = 10_DP**(a6*L+a7 + x234)

    return
  end function ne_ps


  function ne_trough(L, amlt, akp)
    ! (Basically ne_eq_trough, in ne_inner_ps_trough.for)

      ! ---- Trough profile --------------
    !  assuming constant density increase of 0.56+-0.08 cm-3h-1 before phitp
    !  and constant density decrease of -0.83+-0.15 cm-3h-1 after that and before
    !  1 hour MLT, where the decrease is replaced by zero growth.  Growth begins
    !  again at 3.5 hours MLT.  A minimum density of 0.18 cm-3 is used.
    !  Transitions between growth rates are made over an hour or more, as shown below.

    real(kind=DP) :: L, amlt, akp, ne_trough
    real(kind=DP) :: phitp, antp, damping, damping_time
    real(kind=DP) :: down_time, del, center, diff
    real(kind=DP) :: aminden, width, denmin, dengrow
    real(kind=DP) :: sdel, shift, switch0, switch1, switch2, switch3
    real(kind=DP) :: dendamp, geosync_trough, switch_trough2cap, ne_eq_trough


    
    !  compute MLT where trough density peaks assuming constant Kp
    phitp=0.145*akp**2 - 2.63*akp+21.86

    ! print *, 'phitp:', phitp
    !  compute unsmoothed peak density at phitp
    antp=(phitp-3.5)*0.56
    damping_time=min(26.0-phitp,antp/0.83)
    !  given the loss rate above as a minimum.
    damping= -1.0*antp/damping_time
    down_time=phitp+damping_time
    del=3.5-(down_time-24.0)
    center=3.5-del/2.0
    if(center.lt.0.0) center=24.0+center 
    diff=amlt-center
    if(diff.lt.-12.0) diff=24.0+diff
    if(diff.gt.12.0) diff=diff-24.0

    aminden=0.18  ! Minimum density (cm^-3)
    width=2.0*del
    denmin=aminden+diff**2/(del*width)

  !  must deal with computation differently based on current mlt
    dengrow=0.56*(amlt-3.5)+aminden
    sdel=0.4
    shift=0.5
    switch1=switch(amlt,3.5+shift,sdel)
    switch2=switch(amlt,phitp,0.5_DP)

    if(amlt.lt.8.0_DP) then
      dendamp=antp+damping*(amlt+24.0-phitp)
      switch0=switch(amlt,down_time-24.0-shift,sdel)
      geosync_trough= denmin*switch0*(1.0-switch1) +& 
                      &dendamp*(1.0-switch0) +&
                      &dengrow*switch1*(1.0-switch2)
      ! print *,'lt.8.0:',denmin,switch0,switch1,dendamp,dengrow
    else
      dendamp=antp+damping*(amlt-phitp)
      switch3=switch(amlt,down_time-shift,sdel)
      geosync_trough= denmin*switch3 + dengrow*switch1*(1.0-switch2) +&
                      &dendamp*switch2*(1.0-switch3)
      ! print *,'ge.8.0:',denmin,switch3,dengrow,switch1,switch2,dendamp
    end if

    ! scale density from trough density at geosynchronous orbit
    ! to L-shell of interest "l" using a power law of -4.5 in L
    !     =ne_trough[at L=6.6]*al**(-4.5)/6.6**(-4.5)
    ne_eq_trough=geosync_trough*L**(-4.5_DP)/2.0514092e-4_DP



    ! ! Blend into the cap model for higher latitudes
    ! switch_trough2cap = switch(abs(lam), 20.0_DP, 10.0_DP)
    ! ne_trough = ne_eq_trough*(1.0 -switch_trough2cap) + (switch_trough2cap)*ne_polar_cap

    ! No latitude variation
    ne_trough = ne_eq_trough
    

  end function ne_trough



  function check_crossing(a8, amlt, akp, iyear, doy)

      real(kind=DP) :: a8, amlt, akp
      real(kind=DP) :: check_crossing
      real(kind=DP) :: stepl, zl, a, b, c, diff
      integer(kind=DP) :: icount
      real(kind=DP) :: iyear, doy
      stepl=0.5
      ! print *, 'inside check crossing', a8
      zl= a8
      b = pp_profile(zl,amlt,akp,a8)
      a = ne_ps(zl, a8, b, iyear, doy, akp)
      c = ne_trough(zl, amlt, akp)

      diff=a*b - c

      icount=0
      do while (abs(stepl).gt.0.05)
        if ((diff.lt.0.0).and.(stepl.gt.0.0) .or.&
           &(diff.gt.0.0).and.(stepl.lt.0.0)) stepl=-stepl/2.0
        zl=zl+stepl

        b = pp_profile(zl,amlt,akp,a8)
        a = ne_ps(zl, a8, b, iyear, doy, akp)
        c = ne_trough(zl, amlt, akp)

        diff=a*b - c

        icount=icount+1
        if (icount.gt.100) then
          ! temp=pp_profile(zl,amlt,akp,a8)
          print *, 'Failed to find knee in check_crossing'
          stop
        endif
      enddo

      check_crossing=zl 
      return
  end function check_crossing



  ! function ne_iono(lat, mlt, alt)
  !   ! A very simple model of ionosphere electron density.
  !   ! Uses a multiple gaussian curve which was fit to IRI2012 model
  !   ! data at the F2 peak, for dayside (MLT 12) and nightside (MLT 0),
  !   ! for geomagnetic latitudes between -90 and 90 degrees.
  !   ! Uses a sigmoid function to smooth between day and night.
  !   ! Returns electron density in cm^-3

  !   real(kind=DP) :: lat, mlt, alt
  !   ! real(kind=DP), parameter :: coef_nite(11), coef_day(11)
  !   real(kind=DP) :: day, nite
  !   real(kind=DP) :: ne_iono
  !   real(kind=DP), dimension(8) :: coef_nite, coef_day
  !   real(kind=DP), parameter :: mltslope = 0.1_DP ! per-hour rate of change
  !   real(kind=DP) :: s, s1, s2
    

  !   coef_day = (/8.04375253e+05,  -1.98290193e+01,   9.84561228e+00,&
  !               &7.69068420e+05, 1.14717628e+01, 4.07735746e+01, &
  !               &1.85437238e+05,-1.85520036e+03/)
  !   coef_nite =(/-1.32546822e+05, -4.06231821e+01, 2.94414728e+01, &
  !               &1.06125847e+05, 1.13202857e+01,   2.21196102e+01, &
  !               &1.64651016e+05,  -1.57715123e+03/)

  !   day = coef_day(1)*exp(-abs((lat - coef_day(2))/coef_day(3))**2) +&
  !        &coef_day(4)*exp(-abs((lat - coef_day(5))/coef_day(6))**2) +&
  !        &coef_day(7) + coef_day(8)*lat

  !   nite = coef_nite(1)*exp(-abs((lat - coef_nite(2))/coef_nite(3))**2) +&
  !         &coef_nite(4)*exp(-abs((lat - coef_nite(5))/coef_nite(6))**2) +&
  !         &coef_nite(7) + coef_nite(8)*lat

  !   s1 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 18.0_DP)/mltslope))
  !   s2 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 6.0_DP)/mltslope))
  !   s = s1 - s2
  !   ne_iono = s*day + (1.0 - s)*nite
  !   return
  ! end function ne_iono



  function ne_iono(lat, mlt, alt)
  !   ! A very simple model of ionosphere electron density.
  !   ! Uses a multiple gaussian curve which was fit to IRI2012 model
  !   ! at 1500 km, and a density gradient which was calculated as the
  !   ! slope in log space bewteen density at 1500km, and the density
  !   ! at the F2 peak (~350 km +- 100km or so).
  !   ! Uses the fitted gradient to extrapolate densities over altitudes.
  !   ! Uses a sigmoid function to smooth between day and night.
  !   ! Returns electron density in cm^-3.
    real(kind=DP) :: lat, alt, mlt 
    real(kind=DP) :: dens_day, dens_nite, grad_day, grad_nite
    real(kind=DP) :: s, s1, s2, ne_iono
    integer(kind=DP) :: i

    real(kind=DP), parameter :: mltslope = 0.1_DP ! per-hour rate of change
    real(kind=DP), parameter :: dens_coef_day(9) =(/&
                  &9.23183e+03,-2.23382e+01,1.49365e+01,&
                  &1.70763e+04,2.63301e+01,3.75599e+01,&
                  &9.15522e+03,-3.31022e+01,-6.13435e-01/)
    real(kind=DP), parameter :: grad_coef_day(11) = (/&
                  &2.31678e-22,1.17475e-20,-5.20743e-18,&
                  &-1.98686e-16,4.44845e-14,1.19067e-12,&
                  &-1.85079e-10,-3.38517e-09,3.84124e-07,&
                  &3.35202e-06,-1.84164e-03/)
    real(kind=DP), parameter :: dens_coef_nite(9) =(/&
                  &6.99184e+03,-3.11663e+00,1.30464e+01,&
                  &8.58528e+03,2.19513e+01,1.56983e+01,&
                  &2.89385e+03,-1.87291e+01,3.21094e-01/)
    real(kind=DP), parameter :: grad_coef_nite(11) = (/&
                  &2.39859e-23,-1.28908e-20,-1.77123e-20,&
                  &2.32415e-16,-4.42548e-15,-1.46548e-12,&
                  &2.84474e-11,3.59593e-09,-2.62414e-08,&
                  &-2.70750e-06,-1.63765e-03/)

    ! Get densities at 1500 km altitude
    dens_day = dens_coef_day(1)*exp(-abs((lat - dens_coef_day(2))/dens_coef_day(3))**2) +&
              &dens_coef_day(4)*exp(-abs((lat - dens_coef_day(5))/dens_coef_day(6))**2) +&
              &dens_coef_day(7) + dens_coef_day(8)*lat + dens_coef_day(9)*lat**2

    dens_nite = dens_coef_nite(1)*exp(-abs((lat - dens_coef_nite(2))/dens_coef_nite(3))**2) +&
               &dens_coef_nite(4)*exp(-abs((lat - dens_coef_nite(5))/dens_coef_nite(6))**2) +&
               &dens_coef_nite(7) + dens_coef_nite(8)*lat + dens_coef_nite(9)*lat**2

    ! Get gradients (slope in log space to the F2 peak)
    grad_day  = 0
    grad_nite = 0
    do i=1,11
      grad_day  = grad_day +  grad_coef_day(i)*lat**(11-i)
      grad_nite = grad_nite + grad_coef_nite(i)*lat**(11-i)
    end do

    ! Extrapolate in log space to the desired altitude
    dens_day  = dens_day*10.0_DP**(grad_day*(alt - 1500.0_DP))
    dens_nite = dens_nite*10.0_DP**(grad_nite*(alt - 1500.0_DP))

    ! Fade between day and night

    ! if (lat.gt.85.0) then
    !   print *, lat, mlt
    ! end if
    s1 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 18.0_DP)/mltslope))
    s2 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 6.0_DP)/mltslope))
    s = s1 - s2
    ! s = s1
    ! We're done!
    ne_iono = s*dens_day + (1.0 - s)*dens_nite


    ! print *,lat, mlt, alt, dens_day, grad_day, dens_nite, grad_nite, ne_iono
    return
  end function ne_iono



  function ne_cap(lat, r, mlt, akp)

    real(kind=DP) :: lat, r, mlt, akp
    real(kind=DP) :: ne_cap
    real(kind=DP) :: h, refn, powern
    real(kind=DP) :: ne_iono_source
    h = (r - REkm)
    ne_iono_source = ne_iono(lat, mlt, 350.0_DP)
    refn = log(ne_iono_source) + 16.764
    powern = -2.8618
    ne_cap = exp(powern*log(h) + refn) + 0.001
    ne_cap = min(ne_iono_source, ne_cap) 

    ! ne_cap = exp(-3.09*log(h - iono_peak_altitude) + log(ne_iono_source) + 16.764) + 0.001
    ! ne_cap = min(ne_iono_source, ne_cap) 
    return
  end function ne_cap


  subroutine poleward_edge(amlt, akp, edge_lat, edge_L)
    ! calculates the latitude and L-shell of the
    ! boundary between the polar cap region and the rest of
    ! the model.

    real(kind=DP) :: amlt, akp
    real(kind=DP) :: edge_lat, edge_L

    real(kind=DP) :: bmlt, diffmlt
    real(kind=DP) ::diffkp
    integer(kind=DP) :: imlt, jmlt, ikp, jkp
    real(kind=DP) :: pn1, pn2
    real(kind=DP) :: alatcritn, alcrit
    real(kind=DP) :: tranlow, tranhigh
    ! poleward auroral edge
    !  altrans = the half width in L-shell of over which the transition
    !            takes place between the trough and polar cap models.
    !            The L-shell at which this transition is centered is
    !            obtained from PN, which hold empirical locations for
    !            the polarward edge of the auroral zone for several
    !            values of Kp and magnetic local time.

real(kind=DP), parameter :: PN(72, 10) =reshape(&
  &(/73.90,73.95,74.09,74.29,74.51,74.73,74.93,75.10,75.22,75.29,&
  &75.33,75.32,75.27,75.18,75.07,74.95,74.83,74.72,74.63,74.57,&
  &74.53,74.50,74.51,74.59,74.78,75.10,75.55,76.10,76.70,77.30,&
  &77.85,78.35,78.77,79.12,79.41,79.64,79.82,79.96,80.06,80.13,&
  &80.16,80.17,80.15,80.10,80.02,79.92,79.79,79.63,79.45,79.26,&
  &79.05,78.82,78.57,78.30,78.00,77.66,77.27,76.86,76.43,76.00,&
  &75.62,75.32,75.06,74.83,74.60,74.38,74.15,73.95,73.79,73.65,&
  &73.52,73.39,73.30,73.27,73.33,73.47,73.67,73.90,74.13,74.36,&
  &74.57,74.76,74.93,75.07,75.18,75.26,75.32,75.33,75.32,75.27,&
  &75.19,75.10,74.99,74.88,74.79,74.74,74.78,74.93,75.21,75.58,&
  &76.03,76.50,76.98,77.44,77.86,78.23,78.54,78.77,78.94,79.06,&
  &79.12,79.14,79.12,79.07,78.99,78.89,78.76,78.61,78.43,78.23,&
  &78.02,77.81,77.60,77.40,77.18,76.95,76.71,76.46,76.22,75.99,&
  &75.77,75.55,75.36,75.18,75.02,74.85,74.64,74.39,74.11,73.82,&
  &73.54,73.28,73.04,72.81,72.63,72.53,72.53,72.63,72.81,73.05,&
  &73.30,73.55,73.80,74.05,74.29,74.54,74.79,75.02,75.20,75.31,&
  &75.34,75.31,75.24,75.14,75.02,74.91,74.84,74.86,74.98,75.22,&
  &75.56,75.96,76.37,76.75,77.09,77.38,77.62,77.79,77.89,77.93,&
  &77.92,77.87,77.79,77.70,77.59,77.47,77.34,77.20,77.06,76.92,&
  &76.79,76.66,76.54,76.42,76.32,76.22,76.14,76.06,75.98,75.89,&
  &75.79,75.67,75.53,75.37,75.21,75.06,74.91,74.73,74.48,74.17,&
  &73.82,73.47,73.16,72.94,72.81,72.77,72.81,72.91,73.05,73.23,&
  &73.41,73.60,73.80,73.99,74.16,74.31,74.45,74.58,74.70,74.82,&
  &74.90,74.96,74.99,75.00,75.00,75.00,75.01,75.03,75.07,75.12,&
  &75.21,75.32,75.45,75.61,75.77,75.94,76.10,76.25,76.40,76.54,&
  &76.67,76.78,76.87,76.92,76.94,76.91,76.83,76.72,76.58,76.42,&
  &76.27,76.13,76.01,75.90,75.81,75.72,75.64,75.58,75.54,75.51,&
  &75.48,75.45,75.42,75.37,75.32,75.25,75.19,75.13,75.05,74.91,&
  &74.73,74.49,74.21,73.93,73.68,73.51,73.44,73.46,73.56,73.72,&
  &73.91,74.09,74.25,74.40,74.54,74.67,74.79,74.90,74.99,75.07,&
  &75.12,75.14,75.12,75.07,75.00,74.92,74.85,74.79,74.75,74.71,&
  &74.69,74.70,74.75,74.83,74.94,75.07,75.21,75.34,75.46,75.58,&
  &75.68,75.77,75.85,75.91,75.95,75.98,75.99,75.97,75.93,75.87,&
  &75.79,75.69,75.58,75.46,75.34,75.23,75.14,75.06,75.00,74.96,&
  &74.93,74.92,74.93,74.94,74.96,74.98,74.99,74.99,74.98,74.98,&
  &74.96,74.91,74.80,74.64,74.44,74.24,74.07,73.97,73.95,74.01,&
  &74.11,74.26,74.41,74.55,74.68,74.79,74.90,75.00,75.08,75.16,&
  &75.23,75.27,75.29,75.26,75.19,75.08,74.92,74.73,74.52,74.31,&
  &74.11,73.93,73.79,73.71,73.68,73.71,73.79,73.93,74.10,74.31,&
  &74.52,74.73,74.90,75.04,75.12,75.16,75.17,75.16,75.12,75.08,&
  &75.04,74.99,74.94,74.89,74.84,74.79,74.75,74.72,74.71,74.70,&
  &74.70,74.70,74.70,74.70,74.69,74.68,74.66,74.63,74.60,74.57,&
  &74.55,74.57,74.60,74.61,74.58,74.52,74.42,74.32,74.26,74.28,&
  &74.39,74.55,74.75,74.96,75.14,75.28,75.39,75.48,75.55,75.60,&
  &75.63,75.63,75.63,75.61,75.57,75.50,75.39,75.21,74.96,74.65,&
  &74.29,73.91,73.53,73.17,72.88,72.70,72.64,72.68,72.81,73.00,&
  &73.20,73.40,73.59,73.77,73.92,74.04,74.12,74.17,74.19,74.20,&
  &74.21,74.20,74.20,74.20,74.19,74.18,74.16,74.14,74.11,74.09,&
  &74.06,74.04,74.02,74.01,74.00,74.00,74.00,74.00,74.00,74.00,&
  &74.00,74.01,74.05,74.15,74.29,74.44,74.57,74.66,74.73,74.77,&
  &74.82,74.91,75.04,75.17,75.29,75.39,75.46,75.49,75.49,75.48,&
  &75.45,75.41,75.34,75.26,75.18,75.10,75.02,74.94,74.85,74.71,&
  &74.51,74.23,73.87,73.46,73.01,72.56,72.18,71.90,71.74,71.71,&
  &71.77,71.91,72.08,72.26,72.45,72.63,72.80,72.95,73.08,73.20,&
  &73.32,73.42,73.51,73.59,73.64,73.68,73.70,73.72,73.74,73.76,&
  &73.78,73.79,73.79,73.78,73.77,73.74,73.71,73.67,73.62,73.58,&
  &73.56,73.57,73.60,73.67,73.79,73.97,74.18,74.40,74.59,74.74,&
  &74.86,74.95,75.04,75.18,75.35,75.51,75.65,75.72,75.73,75.68,&
  &75.59,75.47,75.34,75.19,75.01,74.81,74.59,74.36,74.11,73.83,&
  &73.52,73.18,72.80,72.40,71.99,71.59,71.19,70.84,70.53,70.30,&
  &70.14,70.06,70.06,70.15,70.29,70.51,70.77,71.06,71.38,71.70,&
  &72.02,72.33,72.61,72.85,73.07,73.26,73.43,73.55,73.63,73.63,&
  &73.56,73.41,73.23,73.02,72.81,72.62,72.45,72.32,72.23,72.18,&
  &72.20,72.28,72.40,72.56,72.76,73.01,73.30,73.66,74.05,74.43,&
  &74.77,75.06,75.29,75.49,75.65,75.79,75.89,75.93,75.92,75.87,&
  &75.79,75.70,75.59,75.47,75.34,75.19,75.01,74.81,74.59,74.36,&
  &74.11,73.83,73.52,73.18,72.80,72.40,71.99,71.59,71.19,70.84,&
  &70.53,70.30,70.14,70.06,70.06,70.15,70.29,70.51,70.77,71.06,&
  &71.38,71.70,72.02,72.33,72.61,72.85,73.07,73.26,73.43,73.55,&
  &73.63,73.63,73.56,73.41,73.23,73.02,72.81,72.62,72.45,72.32,&
  &72.23,72.18,72.20,72.28,72.40,72.56,72.76,73.01,73.30,73.66,&
  &74.05,74.43,74.77,75.06,75.29,75.49,75.65,75.79,75.90,75.95/),(/72,10/))

      oldmlt=amlt
      oldkp=akp
      bmlt=amlt*3.0+1.0
      imlt=int(bmlt)
      diffmlt=bmlt-real(imlt, kind=DP)
      if(imlt.gt.72) imlt=1
      jmlt=imlt+1
      if(jmlt.gt.72) jmlt=1
  !
      ikp=akp+1.0
      diffkp=akp-aint(akp)
      if(ikp.gt.10) ikp=10
      jkp=ikp+1
      if(jkp.gt.10) jkp=10
  !       
      pn1=(pn(jmlt,ikp)-pn(imlt,ikp))*diffmlt+pn(imlt,ikp)
      pn2=(pn(jmlt,jkp)-pn(imlt,jkp))*diffmlt+pn(imlt,jkp)
      !ps1=(ps(jmlt,ikp)-ps(imlt,ikp))*diffmlt+ps(imlt,ikp)
      !ps2=(ps(jmlt,jkp)-ps(imlt,jkp))*diffmlt+ps(imlt,jkp)
  !
  !  Here we must determine the L-shell locations over which the transition
  !  will be made from the trough/plasmasphere model to the polar cap model.
      !alatcrits=(ps2-ps1)*diffkp+ps1
      alatcritn=(pn2-pn1)*diffkp+pn1

      alcrit=1.0/cos(alatcritn*D2R)**2

      ! print *,'auroral zone:',alcrit,altrans,tranlow,tranhigh
      edge_lat = alatcritn
      edge_L = alcrit

  
  end subroutine poleward_edge



  ! subroutine iri_ps_bridge(L)
  !   real(kind=DP) :: eqh    ! Equatorial height
  !   real(kind=DP) :: transh ! Altitude of F2 peak? it's 1.05454 R_E constant.
  !   transh = 350.0 
  !   eqh=(L-1.0)*REkm ! altitude in km

  ! end subroutine iri_ps_bridge








  function main_ps_density(L, a8, a9, iyear, doy, akp, amlt, lam, r)
      real(kind=DP), intent(in) :: L, a8, a9, iyear, doy, akp, amlt, lam, r
      real(kind=DP) :: main_ps_density
      real(kind=DP) :: ne_eq_ps, ne_eq_trough, zl, switch_ps2trough
      real(kind=DP) :: cap_edge_lat, cap_edge_L, ne_polar_cap, switch_cap


      ! find the transition point between the two models
      ! (either the plasmapause, or the point at which
      ! the two curves cross)
      zl = check_crossing(a8, amlt, akp, iyear, doy)

      ! ------- Inner plasmasphere profile -------------
      ne_eq_ps = ne_ps(L, a8, a9, iyear, doy, akp)

      ! ------- Trough (outer plasmasphere) profile -------  
      ne_eq_trough = ne_trough(L, amlt, akp)

   
      if (do_trough.eq.1) then
        switch_ps2trough = switch(L, zl, 0.6_DP) ! sets plasmapause slope...
      else
        switch_ps2trough = 0
      end if

      main_ps_density  = ne_eq_ps*(1.0 - switch_ps2trough) + &
                         &switch_ps2trough*ne_eq_trough
                   
      if (do_cap.eq.1) then
      ! ------- Polar Cap profile ---------------------
        ne_polar_cap = ne_cap(lam, r, amlt, akp)
        call poleward_edge(amlt, akp, cap_edge_lat, cap_edge_L)
        switch_cap = switch(L, cap_edge_L, altrans)
        main_ps_density = main_ps_density*(1.0 - switch_cap) + switch_cap*ne_polar_cap
      end if
      return
  end function main_ps_density



  subroutine find_intersection_iono_ps(a8, a9, iyear, doy, akp, lamr, amlt, merge_alt)
      ! Find the intersection altitude between the ionosphere and plasmasphere
      ! models. Assumes the ionosphere is higher density... I haven't tried it
      ! with other situations.

      real(kind=DP), intent(in) :: a8, a9, iyear, doy, akp, lamr, amlt
      real(kind=DP), intent(out) :: merge_alt

      real(kind=DP) :: alt_guess, iono, iono_lam, ps, L_cur, stepl, diff
      integer(kind=DP) :: i, icount

      stepl = 1000.0_DP     ! Initial step size, km
      alt_guess = 2000.0_DP ! Initial altitude, km

      icount=0
      do while ((abs(stepl).gt.100.0_DP).and.(alt_guess.lt.10000.0_DP))
        if ((diff.lt.0.0).and.(stepl.gt.0.0) .or.&
           &(diff.gt.0.0).and.(stepl.lt.0.0)) stepl=-stepl/2.0
        alt_guess = alt_guess + stepl
  
        
        L_cur = (alt_guess + REkm)/(REkm*cos(lamr)**2.0_DP)
        ps = main_ps_density(L_cur, a8, a9, iyear, doy, akp, amlt, lamr*R2D, (alt_guess + REkm))
        iono = ne_iono(lamr*R2D, amlt, alt_guess)

        diff=iono - ps

        icount=icount+1
        if (icount.gt.100) then ! Should only take 5 or 6...
          ! temp=pp_profile(zl,amlt,akp,a8)
          print *, icount, ps, iono, diff, alt_guess, stepl
          print *, 'Failed to find intersection'
        endif
      enddo

      merge_alt = alt_guess
      return
  end subroutine find_intersection_iono_ps






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

    real(kind=DP) :: ce,ch,che,co
    real(kind=DP) :: x(3),x_gsm(3)
    
    ! Tsyganenko parameters
    integer :: year, day, hour, minute, sec
    real(kind=DP) :: parmod(10)
    integer :: iopt
    ! Tsyganenko corrections
    real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
    real(kind=DP) :: W1, W2, W3, W4, W5, W6
    ! Base B field 
    real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE

    type(simpleStateDataP) :: datap

    real(kind=DP) :: r, lam, lamr
    real(kind=DP) :: p(3)
    real(kind=DP) :: amlt, akp, a8, a9

    real(kind=DP) :: g, h, ne_eq_ps
    real(kind=DP) :: L ! L-shell
    real(kind=DP) :: ne_iono_source, ne_eq_iono, iono_lam, iono_merge_lam, diono
    real(kind=DP) :: iono_merge_L, ne_merge_ps, iono_peak_L
    real(kind=DP) :: iyear, doy, doy_factor, x234
    real(kind=DP) :: ne_polar_cap

    real(kind=DP) :: a6, a7, rz12
    real(kind=DP) :: ne_eq_trough, zl, offset, diff
    real(kind=DP) :: cap_edge_lat, cap_edge_L, switch_cap
    real(kind=DP) :: tranlow, tranhigh

    real(kind=DP) :: al, aheight
    real(kind=DP) :: switch_ps2trough, switch_iono2ps
    real(kind=DP) :: ne_merge_trough, switch_iono2psOrTrough, ne_merge_value
    real(kind=DP) :: switch_trough2cap_lat, switch_trough2cap_rad, eqh

    real(kind=DP) iono_merge_altitude, ne_trough_val, iono_ps_merge_alt

    real(kind=DP) :: aHeH, alphaO, alphaHeP, alphaHe


    ! a6 = -0.79
    ! a7 = 5.208

    iopt = 0

    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)
    ! do_cap    = datap%p%do_cap
    ! do_ps     = datap%p%do_ps
    ! do_iono   = datap%p%do_iono
    ! do_trough = datap%p%do_trough


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
    

    ! ---------- DENSITY MODEL -----------------


    akp = datap%p%kp

    ! r,theta,phi <- x,y,z
    ! theta is azimuth
    ! phi is angle from the z axis
    p = cartesian_to_spherical(x)

    ! MLT from longitude
    if( datap%p%fixed_MLT == 1 ) then
      amlt = datap%p%MLT
    else
      amlt = mod(24.0_DP*p(2)/(2.0_DP*pi)+12.0_DP,24.0_DP)
    end if

    lamr = pi/2.0_DP - p(3)
    lam  = lamr*R2D


    ! ! NOTE lam is in degrees!
    ! lam = max(-89.0_DP, min(89.0_DP, 90.0_DP-(p(3)*R2D)))
    ! lamr = d2r*lam

    ! L-shell (dipole model)
    ! if (R_E*sin(p(3))**2 /= 0.0_DP) then
      L = real(p(1), kind=DP)/(R_E*cos(lamr)**2.0_DP)
  

    ! geocentric radii for all (L,lam) pairs
    r = (REkm) * L * cos(lamr)**2 
  
    ! Location of the plasmapause (Lpp = a8).
    ! a9 is a dropoff rate used in the expression
    ! for H.
    ! print *, 'before pp', a8
    a9=pp_profile(r/REkm,amlt,akp,a8)
    ! print *, 'after pp', a8
    iyear=datap%p%itime(1)/1000
    doy=real(datap%p%itime(1) - iyear*1000, kind=DP)

    ! ------- Polar Cap profile ---------
    ne_polar_cap = ne_cap(lam, r, amlt, akp)

    ! if(oldmlt.ne.amlt .or. oldkp.ne.akp) then
    call poleward_edge(amlt, akp, cap_edge_lat, cap_edge_L)
    ! print *,'mlt: ', amlt, ' cap edge (L): ', cap_edge_L
    tranlow= cap_edge_L - altrans
    tranhigh=cap_edge_L + altrans
    ! end if
    
    ! ! Combined plasmasphere + trough density
    ! ne_eq_ps = main_ps_density(L, a8, a9, iyear, doy, akp, amlt, lamr, r)

    ! find the transition point between the two models
    ! (either the plasmapause, or the point at which
    ! the two curves cross)
    ! print *, 'before check crossiong', a8
    zl = check_crossing(a8, amlt, akp, iyear, doy)
    ! print *, 'zl found is', zl


 

      ! ! ------- Main plasmasphere profile -------------
      ! This includes Carpenter/Anderson PS model, trough model, and cap model
      ne_eq_ps = main_ps_density(L, a8, a9, iyear, doy, akp, amlt, lam, r)
      ! Get altitude where ionosphere and plasmasphere intersect:
      call find_intersection_iono_ps(a8, a9, iyear, doy, akp, lamr, amlt, iono_merge_altitude)
      ! print *, 'after main call', a8
      ne_eq_iono = ne_iono(lamr*R2D, amlt, r - REkm)
      switch_iono2ps = 1.0 - switch(r - REkm, iono_merge_altitude - iono_merge_radius/2.0, iono_merge_radius)
      ne_eq_ps = ne_eq_iono*switch_iono2ps + (1.0 - switch_iono2ps)*ne_eq_ps
    

      ce = ne_eq_ps


    ! -------------- Get plasma constituents from electron density ---------
    ! compute He+ to H+ density ratio in the plasmasphere


      aHeH=10.0**(-1.541 - 0.176*r/REkm + 8.557e-3*f107 -1.458e-5*f107*f107)
    ! Helium concentration drops dramatically with transition to high latitudes
    ! and open field lines
      aHeH=aHeH*(1.0-switch_cap)
    ! compute relative O+ density
      aheight=(r-REkm)
      ! alphaO=0.995/(1.0+(aheight-350.0)**2/281250.0)**3+0.005
      alphaO=0.995/(1.0+(aheight)**2/281250.0)**3+0.005

    ! compute relative He+ concentration in the plasmasphere
      if(aHeH.ne.0.0) then
        alphaHeP=(1.0-alphaO)/(1.0+1.0/aHeH)
        ! alphaHe=max(0.0,alphaHeP*(1.0-exp(-(aheight-400.0)/600.0)))
        alphaHe=max(0.0,alphaHeP*(1.0-exp(-(aheight)/600.0)))
      else
        alphaHe=0.0
      end if

      ce = ce             ! Electrons
      che= alphaHe*ce     ! Helium
      co = alphaO*ce      ! Oxygen
      ch = ce - che - co  ! Hydrogen (whatever's left)
      
    ! end if
    ! species charge and mass
    ! (Electrons, Hydrogen, Helium, Oxygen)
    qs = 1.602e-19_DP*(/ -1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP /);
    ms = (/ 9.10938188e-31_DP, 1.6726e-27_DP, &
         4.0_DP*1.6726e-27_DP, 16.0_DP*1.6726e-27_DP /);
    ! species collision frequencies (currently unused)
    nus = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /);
    Ns = 1.0e6_DP*(/ ce, ch, che, co /);   ! /cm^3  -> /m^3

  
    !  ------------- Tsyganenko magnetic field -------------------
    
    ! Convert the itime parameter into the Tsyganenko parameters
    year = datap%p%itime(1)/1000
    day = mod(datap%p%itime(1),1000)
    hour = datap%p%itime(2)/(1000*60*60)
    minute = (datap%p%itime(2)-hour*(1000*60*60))/(1000*60)
    sec = (datap%p%itime(2)-hour*(1000*60*60)-minute*(1000*60))/(1000)

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
    call tsy_recalc(year, day, hour, minute, sec)
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



end module simple_3d_model_adapter
