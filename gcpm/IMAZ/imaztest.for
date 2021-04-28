c imaztest.for
c
c test program for the imaz subroutine
c
c examples: 
c Geog Lat, Long, Year, DOY, LT, Lv, Ap, h, log(Ne), Uncertainty
c
c Lee-Anne McKinnell, Hermanus Magnetic Observatory, South Africa, April
c 2007
c
      REAL in(9), iout(4,60)

                                                                                
c user input of IMAZ input parameters
c
1       print *,'(geographic)lati/deg,long/deg (lati=100 to end)'
        read(5,*) in(1),in(2)
        if(in(1).gt.90) goto 2
        if(abs(in(1)).lt.60.0.or.abs(in(1)).gt.80.0) then
        	print *,'***Latitude value is outside allowed range',
     &        	 ' (60-80 deg. N or S)'
        	goto 1
        	endif
        print *,'year(yyyy),day(ddd),hour(UT)'
        read(5,*) in(3),in(4),in(5)
        print *,'riometer absorption/dB (-1.0 if not available)'
        read(5,*) in(6)
        print *,'magnetic Ap index (3-hourly)'
        read(5,*) in(7)
        print *,'F10.7 Solar Flux value (daily)'
        read(5,*) in(8)
5        print *,'height/km'
       read(5,*) in(9)
        if(in(9).lt.50.0.or.in(9).gt.150.0) then
        	print *,'***Altitude value is outside allowed range',
     &        	 ' (50-150 km)'
        	goto 5
        	endif

        print *,'  lat  long  year  DOY    UT    Lv    Ap  F10.7  Alt',
     &          '  log(Ne) uncertainty'
        print *,' geog [deg]  YYYY       [hour] [dB]   3-h daily  [km]',
     &          '  [cm-3] '

        call imaz(in,iout)
       if(in(9).ne.-1.0) then
         write(*,120) in(1),in(2),in(3),in(4),in(5),in(6),in(7),
     &     in(8),iout(2,2),iout(3,2),iout(4,2)
       else
         do i=2,60
           write(*,120) in(1),in(2),in(3),in(4),in(5),in(6),in(7),
     &       in(8),iout(2,i),iout(3,i),iout(4,i)
         enddo
       endif

120   format(1X,F5.1,1X,F5.1,2X,F5.0,1X,F4.0,1X,F5.2,1X,F5.2,1X,
     &  F5.1,1X,F5.1,1X,F5.1,1X,F7.4,1X,F7.4)

	goto 1

2	continue
	stop
	end
