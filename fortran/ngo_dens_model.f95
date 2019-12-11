module ngo_dens_model
  use util
  use constants, only : R_E

  implicit real(kind=DP) (a-h)
  implicit real(kind=DP) (o-z)

 real(kind=DP) :: xdra(2000),ydra(2000)
 real(kind=DP) :: l,l0(10),lk,latitu,grarad,gf,latmax,latmin
 
 
 integer :: nreson, ndecre, numres
 real(kind=DP) ::  freq, altinj
 real(kind=DP) :: xtime0(2), xtime1(2)
 integer :: nsuppr
 real(kind=DP) :: a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4)
 real(kind=DP) :: ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk
 real(kind=DP) :: def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,h,hmin
 real(kind=DP) :: pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb
 real(kind=DP) :: rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina
 real(kind=DP) :: therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0
 real(kind=DP) :: hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10)
 real(kind=DP) :: rducun(10),rducus(10),durbe,refalt,hmax
 integer :: init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,numray

contains


  subroutine readinput(filename)
    character(len=*) :: filename
  dimension hducls(10), hducln(10), hducus(10), hducun(10)

!      data pi/3.141592653589793/,clight/2.99792501 e05/,
!     1      consf0/80.61638449/,r0 /6370./

      pi = 3.141592653589793
      clight = 2.99792501e05
      consf0 = 80.61638449
      r0 = 6370.

!
!     constants       (consf0 = q**2 / 4*pi*pi*e0*m)
!-----

      open (unit=4, file=filename, status='old')

! FRF UNUSED
!      zero = 0.0
      radgra = 180./pi
      grarad = 1./radgra
      am = 1.67252e-27 / 9.10955854e-31
      read (4,*) intera, numres, nsuppr, spelat
! FRF UNUSED
      spelat = spelat
      if (intera .eq. 0) goto 50
   50 kinit=2
      init=-1

!      call date (dates)
!      call time (times)
!      xtime0=secnds(0.0)

!!$ dates = fdate()
! times = ctime( time() )
!!$ call ctime( time(), times )
!!$ dummmm = dtime(xtime0)

!     read in the satellite coordinates
!
      count=0.0
!!$      write(*,90)
!!$   90 format (' read input')
      latmax=-100.
      latmin=100.
  105 continue
!     print *, 'past 105'
      read (4,*,end=9501) distre, latitu
      dummmm = 0.0
      if (distre .le. -1.) dummmm = 99999.0
      distre = distre/r0
      if (latitu .le. latmin) latmin=latitu
      if (latitu .ge. latmax) latmax=latitu
  130    continue
            if (dummmm .gt. 99998.0) go to 140
               if (count .gt. 0.0) go to 130
                  if (count .lt. 1.0) go to 105
  140 continue
!     if latmax=latmin, then no satellite trajectory is given, disable
!     variable "nsuppr", latmax and latmin play no role.
      if (abs(latmax-latmin) .le. 1.e-5) nsuppr=0
      latmax=latmax+3.
      latmin=latmin-3.

!     input of model data
      read(4,*,end=9501) num,kskip,mode,kount,kducts,ktape,refalt,dsrrng,dsrlat,dsdens
          if (num .eq. 0) go to 9501

! FRF UNUSED
!          kfile=ktape
          read (4,*,end=9501) egfeq,therm,hm,absb,relb
! FRF UNUSED hm
          hm = hm
          read (4,*,end=9501) rbase,ane0,(alpha0(i),i=2,4)
          read (4,*,end=9501) rzero,scbot,rstop,rdiv,hmin

          if (kducts.eq.0) go to 500
             read (4,*,end=9501) lk,expk,ddk,rconsn,scr

             if (kducts.eq.1) go to 500
                do 935 k=2,kducts
                   read (4,*,end=9501) l0(k),def(k),dd(k),rducln(k),&
                        hducln(k),rducun(k),hducun(k),rducls(k),hducls(k),&
                        rducus(k),hducus(k),sidedu(k)
             hl2n(k)=hducln(k)**2
             hl2s(k)=hducls(k)**2
                hu2n(k)=hducun(k)**2
  935           hu2s(k)=hducus(k)**2
  500 read (4,*,end=9501) pstalt,palt1,palt2,platit,pstlat,plat1,plat2,paltit
!    setting ane0 to desired value at dsrrng,dsrlat
      z(2)=(90.00-dsrlat)*grarad
      z(1)=dsrrng*r0
      call dens
      ane0=ane0*dsdens/ani(1)

      if (pstalt .eq. 0.0) go to 520
!        density profiles
         z(2) = (90.00-platit)*grarad
         z(1) = palt1*r0
         pstep=pstalt*r0
         pfinal=palt2*r0
         i=0
  510    call dens
            palt=z(1)/r0
            z(1)=z(1)+pstep

!           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=palt
            ydra(i)=ani(1)
         if (z(1) .le. pfinal) go to 510

  520 if (pstlat .eq. 0.0) go to 60
         z(1)=paltit+r0
         z(2)=(90.00+plat1)*grarad
         pstep=pstlat*grarad
         pfinal=(90.00+plat2)*grarad
         i=0
  530    call dens
            plat=z(2)*radgra-90.00
            z(2)=z(2)+pstep

!           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=plat
            ydra(i)=ani(1)
         if ( z(2) .le. pfinal) go to 530

   60 continue
 9501 continue
      end subroutine

      

!
      subroutine dens
!    density density density density density density

      dimension exnor(4),qi(4),sh(4)
!
      ! just get rid of "uninitialized variable" warnings due to some
      ! goofy logic below (implied if-else)
      delr = 0.0
      arglr = 0.0
      delroh = 0.0

      cosz2 = cos(z(2))
      sinz2 = sin(z(2))
      cosz22 = cosz2*cosz2
      sinz22 = sinz2*sinz2
      if (init .gt. 0) go to 110
         rb7370 = rbase/7370.
         sh(2) = 1.150600 * therm * rb7370*rb7370
         sh(3) = sh(2)/4.
         sh(4) = sh(3)/4.
         do 105 i = 2,4
  105    alpha(i) = 0.0
!           electron and ion densities.
  110 gph = rbase*(1.0-rbase/z(1))
      exnor(2) = exp(-gph/sh(2))
      exnor(3) = exnor(2)*exnor(2)*exnor(2)*exnor(2)
      exnor(4) = exnor(3)*exnor(3)*exnor(3)*exnor(3)
      q = 0.0
      sumi = 0.0
      do 210 i = 2,num
         qi(i) = alpha0(i)*exnor(i)
         q = q+qi(i)
!        sum1 = -dqdz
         sumi = sumi+qi(i)/sh(i)
  210 continue
      do 211 i = 2,num
  211 alpha(i) = qi(i)/q
!     anr = ne / ne0      radial dependence only - de model
      anr = sqrt(q)
!     anli takes into account the lower ionosphere
      arg = (z(1)- rzero)/scbot
      if (arg.lt.13.0) go to 222
         arg = 13.0
  222 exarg = exp(-arg*arg)
      anli = 1.0-exarg
      dlnlid  = arg*exarg*2./(scbot*anli)
      l = z(1)/(r0*sinz22)
      dlnldr = 0
      dlnldt = 0
      ani(1) = ane0*anr
      if (init .ne. 0) go to 233
         lines = lines+1
  233 ani(1) = ani(1)*anli
      if (kducts.eq.0)go to 250
!        plasmapause
         cotz2 = cosz2/sinz2
         deltal = l-lk
         if (deltal.lt.0.0) go to 410
            d2 = ddk*ddk
            argl = deltal*deltal/(d2*2.0)
            if (argl.lt.80.00) go to 400
               argl = 80.00
  400       f = exp(-argl)
            trm = (rconsn/z(1))**expk
            argr = (z(1)-rconsn)/scr
            if (argr.lt.12.50) go to 405
               argr = 12.5
  405       fr = exp(-argr*argr)
            trmodl = trm+(1.-trm)*fr
            dtrmdr = -expk*trm*(1.-fr)/z(1)-(1.-trm)*fr*2.*argr/scr
            anlk = f+trmodl*(1.0-f)
            factor = deltal*f*l*(1.0-trmodl)/(d2*anlk)
            dlnldr = dlnldr-factor/z(1) + (1.-f)*dtrmdr/anlk
            dlnldt = dlnldt+2.0*factor*cotz2
            ani(1) = ani(1)*anlk
  410    if (kducts.eq.1) go to 250
            if (l0(2) .gt. 0.0) go to 1000
!              sinusoidal density perturbation
               kinit = 3
               deltal = l+l0(2)
               if (deltal*sidedu(2) .ge. 0.0) go to 920
                  deltal = 0
  920          delk = -l0(2)-(lk+ddk)+dd(2)/2
               critl = (lk+ddk)+mod(delk,dd(2))
!              *critl* is the first zero in the sinu pert beyond
!              the knee --  sinusoidal perturbation begins here
               if (l .le. critl) go to 990
                  argl = 2.0*pi*deltal/dd(2)
                  delnl = (def(2)/2.)*(1.+cos(argl))
                  if (latitu.le.0 .and. z(1).le.rducus(2)) go to 950
                  if (latitu.ge.0 .and. z(1).le.rducun(2)) go to 950
!                    continue only if point lies above rducun and rducus
                     if (latitu.ge.0) delr = z(1)-rducun(2)
                     if (latitu.le.0) delr = z(1)-rducus(2)
                     if (latitu.le.0) arglr = delr*delr/hu2s(2)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(2)
                     if (arglr.ge.75.0) go to  990
                        frduct = exp(-arglr)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(2)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(2)
                        go to 960
  950             if (latitu.le.0 .and. z(1).ge.rducls(2)) go to 970
                  if (latitu.ge.0 .and. z(1).ge.rducln(2)) go to 970
!                     continue here only if point lies below 
!                     rducln and rducls
                      if (latitu.ge.0) delr = z(1)-rducln(2)
                      if (latitu.le.0) delr = z(1)-rducls(2)
                      if (latitu.le.0) arglr = delr*delr/hl2s(2)  
                      if (latitu.ge.0) arglr = delr*delr/hl2n(2)
                      if (arglr.ge.75.0) go to  990
                         frduct = exp(-arglr)
                         if (latitu.ge.0) delroh = 2.0*delr/hl2n(2)
                         if (latitu.le.0) delroh = 2.0*delr/hl2s(2)
  960                    delnl = delnl*frduct
                         anl = 1.0+delnl
                         fac = pi*def(2)*l*sin(argl)/(dd(2)*anl)
                         onedut = fac*2.0*cotz2
                         onedur = -fac/z(1) - delnl*delroh/anl
                         go to 980
  970              anl = 1.0+delnl
                   fac = pi*def(2)*l*sin(argl)/(dd(2)*anl)
                   onedut = fac*2.0*cotz2
                   onedur = -fac/z(1)
  980              dlnldr = dlnldr+onedur
                   dlnldt = dlnldt+onedut
                   ani(1) = ani(1)*anl
  990              if (kducts.eq.2) goto 250
!              duct(s)
 1000       do 240 kduc = kinit,kducts
               deltal = l-l0(kduc)
               if (deltal*sidedu(kduc) .ge. 0.0) go to 1200
                  deltal = 0
 1200          d2 = dd(kduc)*dd(kduc)
               argl = deltal*deltal/(d2*2.0)
               if (argl.gt. 80.0 ) go to 240
                  delnl = def(kduc)*exp(-argl)
                  if (latitu.ge.0 .and. z(1).le.rducun(kduc)) goto 1500
                  if (latitu.le.0 .and. z(1).le.rducus(kduc)) goto 1500
!                    continue only if point lies above rducs and rducun
                     if (latitu.ge.0) delr = z(1)-rducun(kduc)
                     if (latitu.le.0) delr = z(1)-rducus(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hu2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = exp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(kduc)
                        go to 1600
 1500             if (latitu.ge.0 .and. z(1).ge.rducln(kduc)) go to 235
                  if (latitu.le.0 .and. z(1).ge.rducls(kduc)) go to 235
!                    continue here only if point lies
!                    below rducln and rducls
                     if (latitu.ge.0) delr = z(1)-rducln(kduc)
                     if (latitu.le.0) delr = z(1)-rducls(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hl2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hl2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = exp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hl2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hl2s(kduc)
 1600                   delnl = delnl*frduct
                        anl = 1.0+delnl
                        fac = delnl*deltal*l/(anl*d2)
                        onedut = fac*2.0*cotz2
                        onedur = -fac/z(1) - delnl*delroh/anl
                        go to 1700
  235             anl = 1.0+delnl
                  onedur = -delnl*deltal*l/(anl*d2*z(1))
                  onedut = -onedur*2.0*z(1)*cotz2
 1700             dlnldr = dlnldr+onedur
                  dlnldt = dlnldt+onedut
                  ani(1) = ani(1)*anl
  240       continue
  250 do 270 i = 2,num
         ani(i) = ani(1)*alpha(i)
  270 continue
!     vzs = dzdr
      vzs = (rbase/z(1))*(rbase/z(1))
!     dlnrdr  = d nr(de) dr / nr
      dlnrdr = -sumi*vzs/(2.0*q)
      do 330 i = 2,num
         dlgndr(i) = -dlnrdr-vzs/sh(i)+dlnldr +dlnlid
         dlgndt(i) = dlnldt
  330 continue
      dlgndr(1) = dlnrdr+dlnldr  +dlnlid
      dlgndt(1) = dlnldt
      return
!    end density end density end density end density
      end subroutine
!
end module ngo_dens_model




