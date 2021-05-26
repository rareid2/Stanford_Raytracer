!
!  This subroutine models the plasmapause location and slope using a
!  modified Lorentzian function as the basis.  It returns a factor
!  that is used elsewhere to effect inclusion of the plasmapause
!  transition from the inner plasmasphere to the trough.
!
!  The factor returned varies from a value of 1 well inside a8
!  to a value of 0 well outside a8.
!
!  The rotation of the bulge with variation with Kpmax has been included
!
!  A plasmapause profile is obtained from Carpenter & Anderson [1991] combined
!  with Higel & Wu [1984] and Moldwin et al, [1994].
!  The profile is included in a8.
!
!  An average Kpmax dependence for the plasmapause slope is also included in
!  a9.
	function pp_profile(al,amlt,akp,a8)
!
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)
	real(kind=DP) :: pp_profile,al,amlt,a8,a9,factor
	real(kind=DP) :: centroid,akp_old,akp,amlt_old
	data akp_old/-99.0/,amlt_old/-99.0/
	print *,'into pp_profile:',al,amlt,akp,a8
!
!  Allow for mlt rotation of the buldge with Kpmax
	if((akp.ne.akp_old) .or. (amlt.ne.amlt_old)) 
     &            call bulge(amlt,akp,a8,a9,centroid)
!	print *,'Recalled subroutine bulge'
	akp_old=akp
	amlt_old=amlt
!
!  compute 10**factor in such a way to avoid floating overflow
	factor=amin1(27.75,2.0*(a9-1.0)*log10(al/a8))
	pp_profile=(1.0+10.0**factor)**(-a9/(a9-1.0))

!     print *,'leaving pp_profile:',pp_profile,factor,a8,a9,centroid
!
	return
	end
