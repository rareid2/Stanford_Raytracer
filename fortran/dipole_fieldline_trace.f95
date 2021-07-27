! This model implements a dipole magnetic field of the earth
module dipole_fieldline_trace
  use util
  use constants
  use types
  use bmodel_dipole
contains

! trace b field in dipole model
function bmodel_fieldline_trace( x, phi )
  ! input in SM coordinates
  real(kind=DP) :: bmodel_fieldline_trace ! this is the output field strength
  real(kind=DP) :: x(3), nx(3)
  real(kind=DP) :: B0tmp(3), Bdir(3)
  real(kind=DP) :: ds, Bmag, direction, stopalt, phi
  integer :: counts
  
  ! get the bfield at the first spot
  B0tmp = bmodel_cartesian( x )
  ! step size in meters
  ds = 100000.0_DP
  ! stop altitude in meters
  stopalt = 400.0e3_DP + R_E 
  ! update position
  nx = x
  
  ! check which way to trace
  if (phi .lt. 0) then
     direction = -1
  else 
     direction = 1
  end if 

  ! check if we are under the stopalt already
  ! then travel up
  r_dist = sqrt(nx(1)*nx(1) + nx(2)*nx(2) + nx(3)*nx(3))
  if (r_dist .lt. stopalt) then
     direction = direction * -1
  end if 

  ! repeat until 1e9 counts reached or stop altitude
  do counts = 1, 1.0e9_DP
     ! print *, counts ! seems like getting into the tens of thousands
     Bmag = sqrt(B0tmp(1)*B0tmp(1) + B0tmp(2)*B0tmp(2) + B0tmp(3)*B0tmp(3))
     Bdir = B0tmp/Bmag 

     ! new position
     nx = (direction * Bdir * ds) + nx

     r_dist = sqrt(nx(1)*nx(1) + nx(2)*nx(2) + nx(3)*nx(3))
     ! check stop condition
     if (r_dist .lt. stopalt) then
        ! print *, 'min alt'
        exit        
     end if
     
     ! update b field at new position
     B0tmp = bmodel_cartesian( nx )

     !if (counts .lt. 5) then
     !   print *, 'x', x 
     !   print *, 'B0tmp', B0tmp
     !   print *, 'Bmag', Bmag
     !   print *, 'nx', nx 
     !end if

  end do
  ! update output
  bmodel_fieldline_trace = Bmag

end function bmodel_fieldline_trace

end module dipole_fieldline_trace
