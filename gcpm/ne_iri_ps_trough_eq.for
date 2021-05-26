!
!	This subroutine is responsible for deteriming the total electron
!	density as a function of position within the ionosphere and
!	plasmasphere at the magnetic equator.
!
!     modified by dlg 1/6/2009 to send itime to ne_inner_ps routine
!
	function ne_iri_ps_trough_eq(al,amlt,akp,itime)
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)
!
	real(kind=DP) :: pi,amltrad,fact,amlt_o,akp_o
	parameter (pi=3.1415927,amltrad=pi/12.0)
	parameter (fact=0.013/pi)
	real(kind=DP) :: al,amlt,akp,ne_iri_ps_trough_eq,ne_inner_ps
	real(kind=DP) :: outf(20,100),oarr(50),a8,ne_eq_trough,re
	real(kind=DP) :: aheight,pp_factor, ps_inner,am1
	real(kind=DP) :: b1,transh,alpha,ano
	real(kind=DP) :: rintercept,ps_bridge,off, pp_profile
	real(kind=DP) :: swtch1,swtch2,switchon,swtch3,alatr,along
	real(kind=DP) :: swtch4,swtch5,ne_eq_trough_den,zl,x234
	real(kind=DP) :: transwidth, check_crossing, diff, offset
	integer(kind=SP) :: itime(2),itime1_o,itime2_o
	data re/6371.0/
	data amlt_o/-99.0/akp_o/-99.0/
	data itime1_o/-99/,itime2_o/-99/
	
!     print *,'entering ne_iri_ps_trough_eq:',al,amlt,akp,itime
!
!  We don't do inside the earth; we're just not that kind of model, humph!
	if(al.le.1.0) then
	  ne_iri_ps_trough_eq=0.0
!       print *,'leaving ne_iri_ps_trough_eq early',ne_iri_ps_trough_eq
	  return
	end if
!
!
!  altitude of point of interest
	aheight=(al-1.0)*re
	pp_factor=pp_profile(al,amlt,akp,a8)
!     print *,'from pp_profile:',pp_factor,al,amlt,akp,a8

	ps_inner=ne_inner_ps(al,amlt,itime,am1,b1,x234)*1.0e6
!     print *,'from ne_inner_ps:',ps_inner,al,amlt,am1,b1,x234

	if(amlt.ne.amlt_o .or. akp.ne.akp_o .or.
     &	itime(1).ne.itime1_o .or. itime(2).ne.itime2_o) then
!
!  determine the height power law fit parameters that connect the
!  topside ionosphere to the plasmasphere
	  call iri_ps_eq_bridge(al,amlt,itime,transh,
     &              alpha,ano,am1,b1,x234,rintercept)

	  amlt_o=amlt
	  akp_o=akp
	  itime1_o=itime(1)
	  itime2_o=itime(2)
	end if
!    print *,'values 1:',a8,am1,b1

	  ps_bridge=ano*aheight**(-alpha)
	!  print *,'bridge-ps n:',rintercept,ano,aheight,alpha
        off=0.0
        transwidth=0.02
	  swtch2=switchon(al,rintercept+off,transwidth)
	  swtch3=switchon(al,rintercept-off,transwidth)

	    alatr=0.0
	    along=(amlt-12.0)*amltrad
!	    along=along - (1.0-sign(1.0,(12.0-amlt)))*pi
	    call iri_sm(alatr,along,al,itime,outf,oarr)
	    swtch1=switchon(aheight,transh,5.0_DP)
	! print *,'entering ne_eq_trough:',swtch1,aheight,transh
          ne_eq_trough_den=ne_eq_trough(al,amlt,akp)
    !      print *,'entering check_crossing:',ne_eq_trough_den
          zl=check_crossing(a8,am1,b1,x234,amlt,akp)
          diff=a8-zl
          offset=(0.0166513 
     &            - 0.0450188*diff)*(1.0-switchon(diff,0.3698744_DP,0.05_DP))
!          offset=amax1(0.0,0.02*((zl-a8-0.2)/0.15))
!          print *,'offset:',akp,a8,zl,(a8-zl),offset
	    swtch4=switchon(al,zl+offset,0.3_DP)
	    swtch5=switchon(al,zl-offset,0.3_DP)
!		print *, 'hey',ps_bridge,swtch2,swtch3
	    ne_iri_ps_trough_eq=outf(1,1)*(1.0-swtch1) +
     &  ((ps_bridge*(1.0-swtch2)*swtch1 +
     &	ps_inner*swtch3)*pp_factor)*(1.0-swtch4) +
     &	ne_eq_trough_den*1.0e6*swtch5
!	print *,'** ',swtch1,swtch2,swtch3,swtch4,swtch5
!	print *,al,rintercept,off,aheight,transh,pp_factor
	! print *,'end of trough_eq:',outf(1,1),ps_bridge,ps_inner,ne_eq_trough_den
    !      print *,ne_iri_ps_trough_eq
!

	return
	end
