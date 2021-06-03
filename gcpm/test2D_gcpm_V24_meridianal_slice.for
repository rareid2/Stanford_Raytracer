c
c	This program exercises this subroutine in 2D

c
	real*4 outn(4),den(201,201),pp_profile
	real*4 r,al,alatr,amlt,akp,ne_inner_ps
	integer*4 itime(2)
	character name*40

      itime(1)=2002185
      ihr=12
      imin=0
      isec=0
      itime(2)=(ihr*3600 + imin*60 + isec) * 1000
	akp=1.0

      do zmlt=0.0,23.0
        mlt_n=zmlt
        mlt_p=zmlt+12.0
        if (mlt_p.ge.24.0) mlt_p=mlt_p-24.0
	do i=0,200
	  x=float(i-100)/10.0
	  amlt=mlt_n
	  if (x .ge. 0.0) amlt=mlt_p
	  type *,'x,zmlt: ',x,zmlt
	do j=0,200
	  z=float(j-100)/10.0
	  r=sqrt(x*x+z*z)
	  den(i+1,j+1)=0.0
	  if (r.gt.1.0) then
	    alatr=atan2(z,abs(x))
          call gcpm_v24(itime,r,amlt,alatr,akp,outn)
          den(i+1,j+1)=outn(1)
        endif
	enddo
      enddo
     
      kp_o=akp
      kp_t=(akp - float(kp_o))*10.0
      write(name,10) mlt_n,mlt_p,kp_o,kp_t
      type *,'writing file: ',name
10    format('gcpm_v24_meridian_',i2.2,'h_',i2.2,'h_kp',
     &       i1.1,'p',i1.1,'.bin')
	open(unit=3,file=name,
     &          type='replace',form='binary')
      write(3) den
      close(unit=3)
      
      enddo

	stop
	end
  