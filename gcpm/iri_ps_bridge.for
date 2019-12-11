c
c	Program runs the subroutine iri_sm to obtain IRI13 densities
c	along an L-shell.
c
c     dlg June 3, 2009  fixed issue with trying to calculate bridge for locations
c                       below the F2 peak along the selected L-shell
c     dlg June 11, 2009 added switchon feature to field aligned bridge function
c                       so that the equatorial density would be reached without
c                       having to use hugh power-law function factors when the
c                       topside fitted power-law was above equatorial density
c                       at the equator.
c
	subroutine iri_ps_bridge(rr,al,alatr,amlt,itime,eq_iri_ps_trough,
     &		transh,rf2,alpha,dno,co,switchh,switchw,istat)
c
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)

	real(kind=DP) :: re,tot_delh,delh,rr,rstart,amltrad
	real(kind=DP) :: rsample1,rsample2
	real(kind=DP) :: r
	parameter (re=6371.0,tot_delh=600.0/re,delh=5.0)
	parameter (r=260.0/re+1.0)
c	parameter (r=260.0/re+1.0,delr=delh/re)
	parameter (amltrad=3.1415927/12.0)
c
	real(kind=DP) :: outf(20,100),oarr(50),alatr,along
	real(kind=DP) :: hs,rf3
	real(kind=DP) :: refden
	real(kind=DP) :: amlt,rs,al
	real(kind=DP) :: ro,transh,alpha,term1
	real(kind=DP) :: eqh,ano,eq_iri_ps_trough,co,fract
	real(kind=DP) :: diffr
	real(kind=DP) :: rloc,rpos,rf2
	real(kind=DP) :: diffold,rlocold,dl,ahemisphere
	real(kind=DP) :: switchh,switchw
	real(kind=DP) :: dno

	real(kind=SP) :: amod
	real(kind=DP) :: cosrl, alatrl, r2, diffh, ah1, ah2, r1, antransh
	real(kind=DP) :: an1, an2, an3

	integer(kind=SP) :: itime(2)
	integer(kind=SP)istat,icount,iflag

	real(kind=DP) :: rz12, f107, neiri, nhoiri, nheiri, noiri
	common /irioutput/ rz12,f107,neiri,nhoiri,nheiri,noiri

d	print *,'entering iri_ps_bridge',rr,al,amlt,itime,eq_iri_ps_trough
c istat must be either 0 or -1 as it is used in an equation later
	istat=0
	dl=al
c !Trevor Garner found error assuming north only, now pass latitude
	ahemisphere=sign(1.0_DP,alatr)
c  get height and densiy of the f2 peak
c	cosrl=amin1(sqrt(r/al),1.0)
c	alatr=acos(cosrl)  !Trevor Garner found error assuming north only, now pass lat
	along=mod(amlt+12.0,24.0_DP)*amltrad
	cosrl=min1(sqrt(r/al),1.0_DP)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r,itime,outf,oarr)
      r2=oarr(2)/re+1.0
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r2,itime,outf,oarr)
      r2=oarr(2)/re+1.0
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r2,itime,outf,oarr)

c approximate the F2 peak along the L-shell=al
	rf2=oarr(2)/re+1.0

c  If L-shell is at or below "r", the starting radial distance
c  for searching for the maximum negative slope above the f2 peak,
c  then the L-shell provided is exclusively an ionospheric issue
c  and we need to pass back parameters that will minimize the hassle
c  associated with the rest of the calculation for density, which
c  will necessarily exclude the bridge density anyway.
d     print *,'f2 peak at:',rf2,al

	if(rr.le.rf2) then
	  istat=-1
d	  print *,'No bridge required, istat=-1 ',rs,al
	  return
	endif

c In an effort to reduce the cals to iri2007 the following is used
c to approximate the point of maximum negative slope in the topside
c ionosphere. This has been obtained from a linear fit to this location
c (derived from the search algorithm above) as a function of returned
c rz12 value from IRI2007. That analysis obtained this relationship:
c     ro = (1.05454+-0.000102) + (8.62678e-5+-1.20975e-6)*rz12
      ro = 1.05454 + 8.62678e-5*rz12
d     print *,'fieldaligned_bridge:',rz12,ro,rf2      
      if (ro .le. rf2) ro=rf2+0.01

	transh=(ro-1.0)*re
      diffh=1.0
      diffr=diffh/re
      ah1=transh-diffh
      ah2=transh+diffh
	r1=ah1/re+1.0
	r2=ah2/re+1.0
c get the density at the maximum slope height
	cosrl=amin1(sqrt(ro/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,ro,itime,outf,oarr)
	antransh=outf(1,1)

c setup for use of densities and heights of the locations
c on either side of the point of maximum negative slope.
c Since only calculating ro from a fitted function, need to separately
c determine the ionospheric densities above and below to support initial
c calculation of the power law function.
	cosrl=amin1(sqrt(r1/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	call iri_sm(alatrl,along,r1,itime,outf,oarr)
	an1=outf(1,1)
d     print *,'an1: ',alatrl,along,r1,an1,al
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	call iri_sm(alatrl,along,r2,itime,outf,oarr)
	an2=outf(1,1)
d     print *,'an2: ',alatrl,along,r2,an2,al

	if(al.le.r2) then
	  istat=-1
d	  print *,'No bridge required, istat=-1 ',al,r2
	  return
	endif

	eqh=(al-1.0)*re
d      print *,'bridge=',ah1,ah2,eqh,transh,antransh
d      print *,'      =',an1,an2,eq_iri_ps_trough
      alpha=-log10(an1/an2)/log10(ah1/ah2)
      ano=an1*ah1**alpha
d      print *,'intial alpha,ano:',alpha,ano
      
      an3=ano*eqh**(-alpha)
d      print *,'setup:',an3,eq_iri_ps_trough

c set up use of switch term that will not function by default
      switchh=eqh*2.0
      switchw=eqh/10.0

      if (eq_iri_ps_trough .ge. an3) then
        if(an2.le.eq_iri_ps_trough) then
d     print *,'inverse IRI-eq:'
          alpha=log10(antransh/eq_iri_ps_trough)/log10(transh/eqh)
          ano=antransh*transh**alpha
          dno=ano
        else
d      print *,'greater than or equal too'
          co=eq_iri_ps_trough - an3
          alpha=-log10((an1-co)/(an2-co))/log10(ah1/ah2)
          ano=(an1-co)*ah1**alpha
          dno=ano
        endif
      else
d      print *,'less than'
c  keep initial alpha and ano values
c  provide switch values that bring the bridge function to the equatorial
c  density at the equator
        switchh=transh+(eqh-transh)/2.0
        switchw=(eqh-transh)/2.0
        dno=ano
        co=0.0
      endif
        
      
d	print *,'final=',dno,alpha,co
d	print *,'    =',ah1,ah2,eqh,an1,an2,eq_iri_ps_trough

d	print *,dens_old,dens,delh,hs
d	print *,'leaving iri_ps_bridge',alpha,ano,transh,switchh,switchw,co
	return
	end
