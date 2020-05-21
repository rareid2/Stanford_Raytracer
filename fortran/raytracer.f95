module raytracer
use types
use util
use constants
use blas
implicit none

real(kind=DP),parameter :: rk45_c(6) = &
     (/0.0_DP,1.0_DP/4.0_DP,3.0_DP/8.0_DP,12.0_DP/13.0_DP,1.0_DP,1.0_DP/2.0_DP/)

real(kind=DP),parameter :: rk45_a2(1) = (/ 1.0_DP/4.0_DP /)
real(kind=DP),parameter :: rk45_a3(2) = &
     (/ 3.0_DP/32.0_DP, 9.0_DP/32.0_DP /)
real(kind=DP),parameter :: rk45_a4(3) = &
     (/ 1932.0_DP/2197.0_DP, -7200.0_DP/2197.0_DP, 7296.0_DP/2197.0_DP /)
real(kind=DP),parameter :: rk45_a5(4) = &
     (/ 439.0_DP/216.0_DP, -8.0_DP, 3680.0_DP/513.0_DP, -845.0_DP/4104.0_DP /)
real(kind=DP),parameter :: rk45_a6(5) = &
     (/ -8.0_DP/27.0_DP, 2.0_DP, -3544.0_DP/2565.0_DP, 1859.0_DP/4104.0_DP, &
       -11.0_DP/40.0_DP /)

real(kind=DP),parameter :: rk45_b4(6) = &
     (/ 25.0_DP/216.0_DP, 0.0_DP, 1408.0_DP/2565.0_DP, 2197.0_DP/4104.0_DP, &
     -1.0_DP/5.0_DP, 0.0_DP /)
real(kind=DP),parameter :: rk45_b5(6) = &
     (/ 16.0_DP/135.0_DP, 0.0_DP, 6656.0_DP/12825.0_DP, &
       28561.0_DP/56430.0_DP, -9.0_DP/50.0_DP, 2.0_DP/55.0_DP /)

contains

! Evaluate the dispersion relation function F(n,w) given the plasma
! parameters and the wavenormal n (in cartesian coordinates)
!
! n = refractive index vector
! w = frequency
! qs = vector of charges
! Ns = vector of number densities in m^-3
! ms = vector of masses
! nus = vector of collision frequencies
! B0 = the magnetic field (vector)
function dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  real(kind=DP) :: dispersion_relation
  real(kind=DP) :: n(3), w, qs(:), Ns(:), ms(:), nus(:), B0(3)
  ! sin^2(phi) and cos^2(phi) respectively
  real(kind=DP) :: sin2phi, cos2phi
  real(kind=DP) :: S,D,P,R,L
  real(kind=DP) :: nmag2, A, B

  ! Find the needed spherical components
  nmag2 = dot_product(n, n)
  ! cos^2(phi)
  cos2phi = (dot_product(n, B0)*dot_product(n, B0)) / &
            (dot_product(n, n)*dot_product(B0, B0))
  sin2phi = 1.0_DP - cos2phi

  ! Find the stix parameters
  call stix_parameters(w, qs, Ns, ms, nus, sqrt(dot_product(B0,B0)), S,D,P,R,L)

  ! Old code
  A = S*sin2phi + P*cos2phi
  B = R*L*sin2phi + P*S*(1.0_DP+cos2phi)
  
  ! If we're sufficiently above the plasma frequency, then treat it 
  ! as free space
  if(w>100.0_DP*sqrt((maxval(Ns)*maxval(abs(qs))**2))/(minval(ms)*EPS0)) then
     dispersion_relation = -nmag2 + 1.0_DP
  else
     dispersion_relation = A*nmag2**2 - B*nmag2 + R*L*P
  end if

!print *, 'A=', A, ', B=', B, ', C=', R*L*P, ', dispers=', dispersion_relation
end function dispersion_relation

! Compute the stix parameters for a multicomponent plasma
! w = frequency
! qs = vector of charges
! Ns = vector of number densities in m^-3
! ms = vector of masses
! nus = vector of collision frequencies
! B0mag = the magnetic field (vector)
subroutine stix_parameters(w, qs, Ns, ms, nus, B0mag, &
                         S,D,P,R,L)
  real(kind=DP), intent(in) :: w, qs(:), Ns(:), ms(:), nus(:), B0mag
  real(kind=DP), intent(out) :: S,D,P,R,L
  real(kind=DP) :: wps2(size(qs)), wcs(size(qs))

  ! Complex raytracing not implemented.  Don't use collisions yet
  ! collisional version
  !!$  wps2 = (Ns*qs**2/ms/EPS0)*(w/(w+j*nus))
  !!$  wcs = ((qs*B0mag)/ms)*(w/(w+j*nus))
  ! noncollisional version
  wps2 = (Ns*qs**2/ms/EPS0)
  wcs = ((qs*B0mag)/ms)

  ! Evaluate the stix parameters given the multicomponent plasma relations
  R = 1.0_DP-sum(wps2/(w*(w+wcs)))
  L = 1.0_DP-sum(wps2/(w*(w-wcs)))
  P = 1.0_DP-sum(wps2/w**2)
  S = 1.0_DP/2.0_DP*(R+L)
  D = 1.0_DP/2.0_DP*(R-L)
end subroutine stix_parameters

! Evaluate the gradient of the dispersion relation with respect to k.
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdk(k, w, x, del, funcPlasmaParams, &
                                  funcPlasmaParamsData)
  real(kind=DP) :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)
  character :: funcPlasmaParamsData(:)
  real(kind=DP) :: dispersion_relation_dFdk(3)
  real(kind=DP) :: d
  real(kind=DP) :: dkx(3), dky(3), dkz(3)

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)

  d = max(del*abs(k(1)), del )
  dkx = d*(/ 1.0_DP,0.0_DP,0.0_DP /)
  dispersion_relation_dFdk(1) = &
       ( dispersion_relation((k+dkx)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dkx)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_DP
  d = max(del*abs(k(2)), del )
  dky = d*(/ 0.0_DP,1.0_DP,0.0_DP /)
  dispersion_relation_dFdk(2) = &
       ( dispersion_relation((k+dky)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dky)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_DP
  d = max(del*abs(k(3)), del )
  dkz = d*(/ 0.0_DP,0.0_DP,1.0_DP /)
  dispersion_relation_dFdk(3) = &
       ( dispersion_relation((k+dkz)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dkz)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_DP

end function dispersion_relation_dFdk


! Evaluate the gradient of the dispersion relation with respect to w
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdw(k, w, x, del, funcPlasmaParams, &
                                  funcPlasmaParamsData)
  real(kind=DP) :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)
  character :: funcPlasmaParamsData(:)
  real(kind=DP) :: dispersion_relation_dFdw
  real(kind=DP) :: d

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)

  d=max( del*abs(w), del )
  
  dispersion_relation_dFdw = &
       ( dispersion_relation(k*c/(w+d), (w+d), qs, Ns, ms, nus, B0 ) - &
         dispersion_relation(k*c/(w-d), (w-d), qs, Ns, ms, nus, B0 ) )/d/2.0_DP

end function dispersion_relation_dFdw


! Evaluate the gradient of the dispersion relation with respect to x.
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdx(k, w, x, del, funcPlasmaParams, & 
                                  funcPlasmaParamsData)
  real(kind=DP) :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)
  character :: funcPlasmaParamsData(:)
  real(kind=DP) :: dispersion_relation_dFdx(3)
  real(kind=DP) :: d
  real(kind=DP) :: dx(3), dy(3), dz(3)
  real(kind=DP) :: n(3)
  real(kind=DP) :: Fp, Fn

  n = k*C/w
  ! Central differencing
  ! x component
  d = max(del*abs(x(1)), del)
  dx = d*(/ 1.0_DP,0.0_DP,0.0_DP /)
  call funcPlasmaParams(x+dx, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dx, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(1) = (Fp-Fn)/d/2.0_DP
  ! y component
  d = max(del*abs(x(2)), del)
  dy = d*(/ 0.0_DP,1.0_DP,0.0_DP /)
  call funcPlasmaParams(x+dy, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dy, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(2) = (Fp-Fn)/d/2.0_DP
  ! dispersion_relation_dFdx(2) = 0.0_DP
  ! z component
  d = max(del*abs(x(3)), del)
  dz = d*(/ 0.0_DP,0.0_DP,1.0_DP /)
  call funcPlasmaParams(x+dz, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dz, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(3) = (Fp-Fn)/d/2.0_DP
end function dispersion_relation_dFdx

! Format of args:
! args(1) = x (meters)
! args(2) = y (meters)
! args(3) = z (meters)
! args(4) = kx (m^-1)
! args(5) = ky (m^-1)
! args(6) = kz (m^-1)
! args(7) = w (rad/s)
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
!
function raytracer_evalrhs(t, args, del, funcPlasmaParams, &
                           funcPlasmaParamsData)
  real(kind=DP) :: args(7), t, del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real(kind=DP) :: x(3), k(3), w, dfdk(3), dfdw, dfdx(3)
  real(kind=DP) :: raytracer_evalrhs(7)
  character :: funcPlasmaParamsData(:)
  
  x = args(1:3)
  k = args(4:6)
  w = args(7)
  
  ! Make del hardcoded for dfdk and dfdw.  dfdx is the one that's really
  ! strongly model-dependent, since some models use real precision.
  dfdk = dispersion_relation_dFdk(k, w, x, 1.0e-8_DP, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdw = dispersion_relation_dFdw(k, w, x, 1.0e-8_DP, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdx = dispersion_relation_dFdx(k, w, x, del, &
                                  funcPlasmaParams, funcPlasmaParamsData)

  raytracer_evalrhs(1:3) = -(dfdk/dfdw)
  raytracer_evalrhs(4:6) = dfdx/dfdw
  raytracer_evalrhs(7) = 0.0_DP
end function raytracer_evalrhs  


! function [stop]=raytracer_stopconditions(pos, k, w, vprel, vgrel, dt)
!
! funcStopConditions should return an error code ~= 0 when some stopping
! criterion is met.  It will take as input position pos, wavenormal k,
! frequency w, relative (scaled by c) phase velocity vprel, relative (scaled
! by c) group velocity vgrel, and the current timestep dt.
!
function raytracer_stopconditions(pos, k, w, vprel, vgrel, dt, nstep, &
                                  maxsteps, minalt)
  real(kind=DP) :: pos(3), k(3), w, vprel(3), vgrel(3), dt, minalt
  integer :: raytracer_stopconditions, nstep, maxsteps

  !print *,' k=', k, ', w=', w, ', vprel=', vprel, ', vgrel=', vgrel
  raytracer_stopconditions = 0
  if( sqrt(dot_product(pos,pos)) < minalt ) then
    ! Minimum altitude reached
    raytracer_stopconditions = 1
    print *, '  Stopping integration.  Specified minimum altitude reached.'
  elseif( sqrt(dot_product(k,k)) == 0.0_DP ) then
    ! Nonsensical k
    raytracer_stopconditions = 2
    print *, '  Stopping integration.  k=0.'
  elseif( sqrt(dot_product(vgrel,vgrel)) > 1.0_DP+1e-2_DP ) then
    ! Faster than light group velocity, with fudge for roundoff error
    raytracer_stopconditions = 3
    print *, '  Stopping integration.  Nonsensical group velocity = ', &
         sqrt(dot_product(vgrel,vgrel))
  elseif( dt < 1e-12 ) then
  ! elseif( dt < 1e-10 ) then
    ! dt too small
    print *, '  Stopping integration.  dt too small.'
    raytracer_stopconditions = 5
  elseif( nstep >= maxsteps ) then
    ! maxsteps exceeded
    print *, '  Stopping integration.  Maxsteps exceeded.'
    raytracer_stopconditions = 6
  end if
end function raytracer_stopconditions

subroutine form_dispersion_matrix(M, n2, phi, S, D, P)
  complex(kind=DP),intent(out) :: M(3,3)
  real(kind=DP),intent(in) :: S,D,P
  real(kind=DP),intent(in) :: n2
  real(kind=DP),intent(in) :: phi

  M(1,1) = cmplx(S-n2*cos(phi)**2, 0)
  M(1,2) = cmplx(0, -D)
  M(1,3) = cmplx(n2*cos(phi)*sin(phi), 0)
  M(2,1) = cmplx(0, D)
  M(2,2) = cmplx(S-n2, 0)
  M(2,3) = cmplx(0,0)
  M(3,1) = cmplx(n2*cos(phi)*sin(phi), 0)
  M(3,2) = cmplx(0,0)
  M(3,3) = cmplx(P-n2*sin(phi)**2, 0)

end subroutine form_dispersion_matrix

function is_right_handed(n2, phi, S, D, P)
  real(kind=DP),intent(in) :: S,D,P
  real(kind=DP),intent(in) :: n2
  real(kind=DP),intent(in) :: phi
  logical :: is_right_handed
  
  complex(kind=DP) :: M(3,3), U(3,3), VT(3,3)
  real(kind=DP) :: svd(3)
  complex(kind=DP) :: E(3)
  real(kind=DP) :: E0(3), E90(3)
  real(kind=DP) :: angle

  call form_dispersion_matrix(M, n2, phi, S, D, P)
  call csvd(M, U, svd, VT)
  E = VT(3,:)

  E0 = real(E*cmplx(1,0,kind=DP))
  E90 = real(E*cmplx(0,1,kind=DP))
  angle = (atan2(E90(2),E90(1)) - atan2(E0(2),E0(1)))

  if( angle > PI ) then
     angle = angle - 2.0_DP*PI
  else if( angle < -PI ) then
     angle = angle + 2.0_DP*PI 
  end if

  if( angle < 0.0_DP ) then
     is_right_handed = .false.
  else
     is_right_handed = .true.
  end if

end function is_right_handed


subroutine solve_dispersion_relation(k, w, x, k1, k2, &
                                     funcPlasmaParams, funcPlasmaParamsData)
! Solve the dispersion relation and return k.  Inputs:
! 
! k = wavenormal
! w = frequency
! x = position
!
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
! k1 and k2 are the positive and negative roots, respectively.
  real(kind=DP),intent(in) :: k(3), w, x(3)
  complex(kind=DP),intent(out) :: k1, k2
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)
  character :: funcPlasmaParamsData(:)
  complex(kind=DP) :: discriminant, n1, n2, nsquared1, nsquared2
  real(kind=DP) :: cos2phi, sin2phi, B0mag
  real(kind=DP) :: S,D,P,R,L, A,B
  real(kind=DP) :: phi

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  
  ! Find the angle the k vector makes with B0
  ! cos^2(phi)
  cos2phi = (dot_product(k, B0)*dot_product(k, B0)) / &
            (dot_product(k, k)*dot_product(B0, B0))
  sin2phi = 1.0_DP - cos2phi
  phi = acos(sqrt(cos2phi))
  ! print *, 'kvectorriley: ', phi

  ! Magnitude of B0
  B0mag = sqrt(dot_product(B0,B0))
  
  ! Find the stix parameters
  call stix_parameters(w, qs, Ns, ms, nus, B0mag, S,D,P,R,L)
  
  ! (Stix, "Waves In Plasmas", page 9)
  A = S*sin2phi+P*cos2phi
  B = R*L*sin2phi+P*S*(1.0_DP+cos2phi)
  discriminant = B**2-4.0_DP*A*R*L*P
  nsquared1 = (B+sqrt(discriminant))/(2.0_DP*A)
  nsquared2 = (B-sqrt(discriminant))/(2.0_DP*A)

  ! print *, 'Bomag: ', B0mag
  ! print *, 'cos2phi: ', cos2phi
  ! print *, 'kmag init: ', dot_product(k, k)
  

  ! print *, 'A: ', A
  ! print *, 'B: ', B
  ! print *, 'nsqared1: ', nsquared1
  ! print *, 'nsquared2 ', nsquared2
  n1 = sqrt(nsquared1)
  n2 = sqrt(nsquared2)
  
  ! default case
  k1 = w*n1/C
  k2 = w*n2/C

  ! Sort modes based on handedness, but only consider handedness of propagating modes
  ! Convention: k2 is the right-handed mode, k1 is the left
  ! print *, 'n1: ', n1
  ! print *, 'n2: ', n2

  if( real(n1) > 0.0_DP .and. is_right_handed(real(nsquared1), phi, S, D, P) ) then
     k1 = w*n2/C
     k2 = w*n1/C
  ! else if( real(n2) > 0.0_DP .and. is_right_handed(real(nsquared2), phi, S, D, P) ) then
  !    k1 = w*n1/C
  !    k2 = w*n2/C
  end if

  ! if ( (real(n1) == 0.0_DP) .and. (real(n2) == 0.0_DP) ) then
  ! ! Both modes are evanescent. Return zero for k vectors instead of NaNs.
  !   print *, 'Evanescent mode'
  !   k1 = 0.0_DP
  !   k2 = 0.0_DP
  ! end if

end subroutine solve_dispersion_relation

function rk4(t, x, del, dt, funcPlasmaParams, funcPlasmaParamsData)
  real(kind=DP) :: x(:), t, dt, del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  character :: funcPlasmaParamsData(:)

  real(kind=DP) :: k1(size(x)),k2(size(x)),k3(size(x)),k4(size(x))
  real(kind=DP) :: rk4(size(x))

  ! Integrate in time using rk4
  ! x' = f(t,x)
  k1 = dt*raytracer_evalrhs(t,x, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k2 = dt*raytracer_evalrhs(t+1.0_DP/2.0_DP*dt, x+1.0_DP/2.0_DP*k1, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k3 = dt*raytracer_evalrhs(t+1.0_DP/2.0_DP*dt, x+1.0_DP/2.0_DP*k2, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k4 = dt*raytracer_evalrhs(t+dt,x+k3, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  rk4 = x + 1.0_DP/6.0_DP*(k1+2.0_DP*k2+2.0_DP*k3+k4)

end function rk4

subroutine rk45(t, x, del, dt, funcPlasmaParams, funcPlasmaParamsData, &
                out4, out5)
  real(kind=DP),intent(in) :: x(:), t, dt, del
  real(kind=DP),intent(out) :: out4(size(x)), out5(size(x))
  complex(kind=DP) :: k1mag, k2mag, k(3)


  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  character :: funcPlasmaParamsData(:)

  real(kind=DP) :: k1(size(x)),k2(size(x)),k3(size(x)),k4(size(x))
  real(kind=DP) :: k5(size(x)),k6(size(x))
  real(kind=DP) :: tmpx(size(x))

  ! Integrate in time using rk45 (embedded 4,5 scheme)
  ! x' = f(t,x)
  tmpx = x
  k1 = dt*raytracer_evalrhs(t,tmpx, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  tmpx = x + ( rk45_a2(1)*k1 )
  k2 = dt*raytracer_evalrhs(&
         t + rk45_c(2)*dt, &
         tmpx, &
         del, funcPlasmaParams, funcPlasmaParamsData)
  tmpx = x + ( rk45_a3(1)*k1 + rk45_a3(2)*k2 )
  k3 = dt*raytracer_evalrhs(&
         t + rk45_c(3)*dt, &
         tmpx, &
         del, funcPlasmaParams, funcPlasmaParamsData)
  tmpx = x + ( rk45_a4(1)*k1 + rk45_a4(2)*k2 + rk45_a4(3)*k3 )
  k4 = dt*raytracer_evalrhs(&
         t + rk45_c(4)*dt, &
         tmpx, &
         del, funcPlasmaParams, funcPlasmaParamsData)
  tmpx = x + ( rk45_a5(1)*k1 + rk45_a5(2)*k2 + rk45_a5(3)*k3 + &
               rk45_a5(4)*k4 )
  k5 = dt*raytracer_evalrhs(&
         t + rk45_c(5)*dt, &
         tmpx, &
         del, funcPlasmaParams, funcPlasmaParamsData)
  tmpx = x + ( rk45_a6(1)*k1 + rk45_a6(2)*k2 + rk45_a6(3)*k3 + &
               rk45_a6(4)*k4 + rk45_a6(5)*k5 )
  k6 = dt*raytracer_evalrhs(&
         t + rk45_c(6)*dt, &
         tmpx, &
         del, funcPlasmaParams, funcPlasmaParamsData)

  out4 = x + &
       (rk45_b4(1)*k1 + rk45_b4(2)*k2 + rk45_b4(3)*k3 + &
        rk45_b4(4)*k4 + rk45_b4(5)*k5 + rk45_b4(6)*k6)
  out5 = x + &
       (rk45_b5(1)*k1 + rk45_b5(2)*k2 + rk45_b5(3)*k3 + &
        rk45_b5(4)*k4 + rk45_b5(5)*k5 + rk45_b5(6)*k6)

end subroutine rk45


! function norm(v, s)
!   integer(kind=DP) :: s
!   real(kind=DP) :: v(s)
!   real(kind=DP) :: norm

!   norm = sqrt(dot_product(v,v))
!   return
! end function norm


subroutine raytracer_run( pos,time,vprel,vgrel,n,&
     B0, qs, ms, Ns, nus, stopcond, &
     pos0, dir0, w0, dt0, dtmax, maxerr, maxsteps, minalt, &
     root, tmax, fixedstep, &
     del, funcPlasmaParams, funcPlasmaParamsData, funcStopConditions)
  
  real(kind=DP), allocatable, intent(out) :: & 
       pos(:,:), time(:), vprel(:,:), vgrel(:,:), n(:,:), & 
       B0(:,:), qs(:,:), ms(:,:), Ns(:,:), nus(:,:)
  real(kind=DP), allocatable :: tmpsize2(:,:), tmpsize1(:)
  integer, intent(out) :: stopcond
  real(kind=DP), intent(in) :: pos0(3), w0, dt0, dtmax, maxerr, tmax
  real(kind=DP), intent(inout) :: dir0(3)
  integer, intent(in) :: root, fixedstep, maxsteps
  real(kind=DP), intent(in) :: del, minalt

  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       use types
       real(kind=DP) :: x(3)
       real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real(kind=DP) :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  interface
     integer function funcStopConditions( pos, k, w, vprel, vgrel, dt, &
          nstep, maxsteps, minalt)
       use types
       real(kind=DP) :: pos(3), k(3), w, vprel(3), vgrel(3), dt, minalt
       integer :: nstep, maxsteps
     end function funcStopConditions
  end interface
  character :: funcPlasmaParamsData(:)

  real(kind=DP) :: dt, x(7), t, x0(7), dfdk(3), dfdw, est1(7), est2(7)
  real(kind=DP) :: dfdk_est1(3), dfdk_est2(3)
  ! Temporary variables needed by call to get and output the plasma
  ! parameters at our next calculated point.
  real(kind=DP), allocatable :: qstmp(:), Nstmp(:), mstmp(:), nustmp(:)
  real(kind=DP) :: B0tmp(3), p(3)

  integer :: lastrefinedown, nstep
  complex(kind=DP) :: k1mag, k2mag, k0mag, k(3), k0(3)
  real(kind=DP) :: dtincr, err, cur_pos(3), w
  real(kind=DP) :: er1, er2, B0mag

  real(kind=DP) :: kmag, kmag_prev, k_err
  real(kind=DP) :: cur_dir(3), prev_dir(3)
  real(kind=DP) :: angle_error, n_angle_refines


  ! Convention: If dir0 is all zeroes, start aligned with the field:
  if ((dir0(1).eq.0).and.(dir0(2).eq.0).and.(dir0(3).eq.0)) then
    call funcPlasmaParams(pos0, qstmp, Nstmp, mstmp, nustmp, B0tmp, funcPlasmaParamsData)
    p = cartesian_to_spherical(pos0)  ! Spherical starting direction

    B0tmp = cartesian_to_spherical_vec(B0tmp, p(2), p(3))
    print *, 'sph_p: ',p
    print *, 'sph_B: ',B0tmp
    b0tmp(1) = abs(b0tmp(1))
    B0tmp = spherical_to_cartesian_vec(B0tmp, p(2), p(3))

    dir0 = B0tmp/sqrt(dot_product(B0tmp,B0tmp))

    print *,'Field-aligned starting direction: ', dir0
  end if







  ! Find k at the given direction
  call solve_dispersion_relation( dir0, w0, pos0, k1mag, k2mag, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  if( root == 1 ) then
     k0mag = k1mag
  else
     k0mag = k2mag
  end if
  k0 = k0mag*dir0

  ! Our state vector -- position, k, and w respectively
  x0 = (/ pos0, real(k0,kind=DP), w0 /)

  dt = dt0
  x = x0
  t = 0.0_DP
  lastrefinedown = 0

  dfdk = dispersion_relation_dFdk(x(4:6), w0, x(1:3), 1.0e-8_DP, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdw = dispersion_relation_dFdw(x(4:6), w0, x(1:3), 1.0e-8_DP, &
                                  funcPlasmaParams, funcPlasmaParamsData)

  ! Find the plasma parameters at our starting point
  call funcPlasmaParams(x(1:3), qstmp, Nstmp, mstmp, nustmp, B0tmp, &
                        funcPlasmaParamsData)

  ! Initialize the output variables
  allocate(pos(3,1))
  pos(:,1) = x(1:3)
  allocate(time(1))
  time(1) = t
  allocate(n(3,1))
  n(:,1) = x(4:6)*C/w0
  allocate(vprel(3,1))
  ! vprel(:,1) = n(:,1)/dot_product(n(:,1),n(:,1))
  allocate(vgrel(3,1))
  ! vgrel(:,1) = -(dfdk/dfdw)/C

  ! If the ray starts off in a non-propagating mode,
  ! check the magnitude of n and avoid writing NaNs
  ! in the file. Stopconds will force an exit on the first
  ! timestep anyway.
  if (dot_product(n(:,1), n(:,1)) > 0) then
    vprel(:,1) = n(:,1)/dot_product(n(:,1),n(:,1))
    vgrel(:,1) = -(dfdk/dfdw)/C
  else 
    vprel(:,1) = 0.0_DP
    vgrel(:,1) = 0.0_DP
  end if
  allocate(B0(3,1))
  B0(:,1) = B0tmp
  allocate(qs(size(qstmp),1))
  qs(:,1) = qstmp
  allocate(ms(size(mstmp),1))
  ms(:,1) = mstmp
  allocate(Ns(size(Nstmp),1))
  Ns(:,1) = Nstmp
  allocate(nus(size(nustmp),1))
  nus(:,1) = nustmp
  stopcond = 0

  print *,'Starting position: ',x(1:3)
  print *,'Starting kay: ',x(4:6)

  nstep = 1
  do
     if( t >= tmax ) then
        ! Normal exit
        stopcond = 0
        exit
     end if
     ! Check stop conditions
     stopcond = funcStopConditions( x(1:3), x(4:6), x(7), &
                                    vprel(:,size(vprel,2)), &
                                    vgrel(:,size(vgrel,2)), dt, &
                                    nstep, maxsteps, minalt ) 
     if( stopcond /= 0 ) then
        print *, 'stopped at t=', t
        exit
     end if

     ! Print the current time step -- comment this out if you don't want
     ! all that output on your screen
     ! print '(a,f5.2)', 't=', t


     if( fixedstep == 0 ) then
        ! Adaptive timesteps - use embedded rk45 scheme
        call rk45( t, x, del, dt, &
             funcPlasmaParams, funcPlasmaParamsData, est1, est2 )

        dtincr = dt
    
        ! Compute the derivatives with respect to k too
        dfdk_est1 = dispersion_relation_dFdk(est1(4:6), w, est1(1:3), &
             1.0e-8_DP, funcPlasmaParams, funcPlasmaParamsData)
        dfdk_est2 = dispersion_relation_dFdk(est2(4:6), w, est2(1:3), &
             1.0e-8_DP, funcPlasmaParams, funcPlasmaParamsData)

        ! Error term is the max of the relative errors in dfdk AND k
        ! This provides better error control when magnetospherically reflecting
        ! (RKF45 solution -- aps, 12.2016)
        err = max( &
             sum(abs(est1(4:6)-est2(4:6))) / sum(abs(est2(4:6))), &
             sum(abs(dfdk_est1-dfdk_est2))/sum(abs(dfdk_est2)))

        if( err > maxerr ) then
          ! if ( dt > 1e-6) then
           ! retry
           ! (Constraining to a minimum timestep... 2.2017)
           !print *, 'Refine down'
           dt=0.8_DP*dt
           ! Prevent refinement loops
           lastrefinedown = 1
           cycle
          ! else
            ! print *, "dt at minimum"
            ! dt=1e-6
          ! end if
        end if
        if( lastrefinedown==0 .and. &
           ! err < maxerr/10.0_DP .and. &
           err < maxerr/100.0_DP .and. &     ! Trying 100 instead of 10 -- 1.6.2018 -- hopefully the larger error tolerance will still have a generally small error in the healthy regions
           dt*1.25_DP < dtmax ) then
           !print *, 'Refine up'
           ! Refine up
           dt=dt*1.25_DP
           lastrefinedown = 0
        end if
     else
        ! Fixed timesteps
        est2 = rk4(t, x, del, dt, funcPlasmaParams, funcPlasmaParamsData)
        dtincr = dt
     end if

     cur_pos = est2(1:3)
     k = est2(4:6)
     w = est2(7)
     ! Refine both estimates based on the physics (must satisfy dispersion
     ! relation)
     call solve_dispersion_relation(real(k,kind=DP), w, cur_pos,k1mag,k2mag, &
                                    funcPlasmaParams, funcPlasmaParamsData )
     if( root == 1 ) then
        ! Preserve direction
        k = k1mag*(k/sqrt(dot_product(k,k)))
        kmag = k1mag
     else
        ! Preserve direction
        k = k2mag*(k/sqrt(dot_product(k,k)))
        kmag = k2mag
     end if

     ! ! refine timestep if magnitude of K-vector changed too much 
     ! if (t > 0) then
     !    k_err = abs((kmag - kmag_prev)/kmag)
     !    kmag_prev = kmag
     !    print *,'t: ',t,' dk: ',k_err
     !  else
     !    k_err = 0.0_DP
     ! end if

     ! if( k_err > 0.1 ) then
     !    if( fixedstep == 0 ) then
     !       print *, 'Refine down: dk over threshold'
     !       dt=dt/2.0_DP
     !       lastrefinedown = 1
     !       cycle        
     !    end if
     ! end if

!------------------------------------------

    !  ! refine timestep if direction of propagation vector changed too much
    !  if (nstep > 1) then
    !   cur_dir = cur_pos - pos(:,nstep)
    !   prev_dir= pos(:, nstep) - pos(:,nstep-1)
    !   angle_error = (180.0_DP/3.14159_DP)*acos(dot_product(cur_dir, prev_dir)/ &
    !              sqrt(dot_product(cur_dir, cur_dir))/ &
    !              sqrt(dot_product(prev_dir, prev_dir)))
      
    ! else
    !   n_angle_refines = 0
    !   angle_error = 0.0_DP
    ! end if
    ! ! print *,'t: ',t,' dTheta: ',angle_error

    ! if (angle_error > 10) then
    !   if (n_angle_refines < 5) then
    !     print *, 'outside angle bounds - refine down'
    !     dt = dt/2.0_DP
    !     n_angle_refines = n_angle_refines + 1
    !     lastrefinedown = 1
    !     cycle
    !   else
    !     print *, 'Maximum refine attempts -- continuing anyway'
    !     n_angle_refines = 0
    !     continue
    !   end if
    ! end if

!-----------------------------------------

     ! refine timestep if outside resonance cone
     if( dot_product(imag(k),imag(k)) > 0.0_DP ) then
        if( fixedstep == 0 ) then
           ! Force a refinement if we've popped outside the resonance cone.
           !print *, 'Refine down: outside of resonance cone'
           dt=dt/2.0_DP
           lastrefinedown = 1
           cycle
        else
           print *, &
                'Cannot continue with fixed timestep - outside resonance cone.'
           print *, 'Try reducing timestep or using adaptive timestepping.'
           return
        end if
     end if
  
     ! Update
     x = est2
     x(4:6) = real(k,kind=DP)
     lastrefinedown = 0
     t = t+dtincr
     nstep = nstep + 1
     
     ! Group velocity
     dfdk = dispersion_relation_dFdk(x(4:6), w, x(1:3), 1.0e-8_DP, &
          funcPlasmaParams, funcPlasmaParamsData)
     dfdw = dispersion_relation_dFdw(x(4:6), w, x(1:3), 1.0e-8_DP, &
          funcPlasmaParams, funcPlasmaParamsData)
     ! Find plasma parameters at our new point
     call funcPlasmaParams(x(1:3), qstmp, Nstmp, mstmp, nustmp, B0tmp, &
                           funcPlasmaParamsData)

!!$     print *, 'pos=',  x(1:3)
!!$     print *, 'time=',  t
!!$     print *, 'n=',  x(4:6)
!!$     print *, 'pos=', x(1:3)
!!$     print *, 'Ns=', Nstmp
!!$     print *, 'dfdk=', dfdk
!!$     print *, 'dfdw=', dfdw

!!$  fprintf('t=%3.3g, x=(%3.3g,%3.3g,%3.3g), n=(%3.3g,%3.3g,%3.3g), vpr=(%3.3g,%3.3g,%3.3g), vgr=(%3.3g,%3.3g,%3.3g)\n', time(end), pos(1,end), pos(2,end), pos(3,end), n(1,end), n(2,end), n(3,end), vprel(1,end), vprel(2,end), vprel(3,end), vgrel(1,end), vgrel(2,end), vgrel(3,end));


     ! Reallocate and update the output variables
     ! pos = [pos, x(1:3)]
     allocate(tmpsize2(size(pos,1), size(pos,2)+1))
     tmpsize2(1:size(pos,1),1:size(pos,2)) = pos
     call move_alloc(tmpsize2, pos)
     pos(:,size(pos,2)) = x(1:3)
     ! time = [time, t]
     allocate(tmpsize1(size(time)+1))
     tmpsize1(1:size(time)) = time
     call move_alloc(tmpsize1, time)
     time(size(time)) = t
     ! n = [n, (x(4:6)*C/w)]
     allocate(tmpsize2(size(n,1), size(n,2)+1))
     tmpsize2(1:size(n,1),1:size(n,2)) = n
     call move_alloc(tmpsize2, n)
     n(:,size(n,2)) = x(4:6)*C/w
     ! vprel = [vprel,  n(:,end)/(norm(n(:,end))^2)]
     allocate(tmpsize2(size(vprel,1), size(vprel,2)+1))
     tmpsize2(1:size(vprel,1),1:size(vprel,2)) = vprel
     call move_alloc(tmpsize2, vprel)
     vprel(:,size(vprel,2)) = & 
          n(:,size(n,2))/dot_product(n(:,size(n,2)),n(:,size(n,2)))
     ! vgrel = [vgrel, -(dfdk./dfdw)/C]
     allocate(tmpsize2(size(vgrel,1), size(vgrel,2)+1))
     tmpsize2(1:size(vgrel,1),1:size(vgrel,2)) = vgrel
     call move_alloc(tmpsize2, vgrel)
     vgrel(:,size(vgrel,2)) = -(dfdk/dfdw)/C
     ! B0 = [B0, B0tmp]
     allocate(tmpsize2(size(B0,1), size(B0,2)+1))
     tmpsize2(1:size(B0,1),1:size(B0,2)) = B0
     call move_alloc(tmpsize2, B0)
     B0(:,size(B0,2)) = B0tmp
     ! qs = [qs, qstmp]
     allocate(tmpsize2(size(qs,1), size(qs,2)+1))
     tmpsize2(1:size(qs,1),1:size(qs,2)) = qs
     call move_alloc(tmpsize2, qs)
     qs(:,size(qs,2)) = qstmp
     ! ms = [ms, mstmp]
     allocate(tmpsize2(size(ms,1), size(ms,2)+1))
     tmpsize2(1:size(ms,1),1:size(ms,2)) = ms
     call move_alloc(tmpsize2, ms)
     ms(:,size(ms,2)) = mstmp
     ! Ns = [Ns, Nstmp]
     allocate(tmpsize2(size(Ns,1), size(Ns,2)+1))
     tmpsize2(1:size(Ns,1),1:size(Ns,2)) = Ns
     call move_alloc(tmpsize2, Ns)
     Ns(:,size(Ns,2)) = Nstmp
     ! nus = [nus, nustmp]
     allocate(tmpsize2(size(nus,1), size(nus,2)+1))
     tmpsize2(1:size(nus,1),1:size(nus,2)) = nus
     call move_alloc(tmpsize2, nus)
     nus(:,size(nus,2)) = nustmp
     
!!$  print '(es24.15e3, 3es24.15e3, 3es24.15e3, 3es24.15e3, 3es24.15e3)', t, & 
!!$       pos(:,size(pos,2)), n(:,size(n,2)), vprel(:,size(vprel,2)), &
!!$       vgrel(:,size(vgrel,2))
       

  end do;

end subroutine raytracer_run
   


end module raytracer
