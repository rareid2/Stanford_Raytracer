!
!	Function varies from 0 to 1 to within 0.1%
!	as the arguement x passes from (a-da) to (a+da).
!	So "da" is the range over which the function turns on.
!
	function switchon(x,a,da)


	real*8 x,a,da,c,switchon
!
	c=3.4534/da
	switchon=tanh(c*(x-a))/2+0.5
	!print *, 'i am groot', switchon,c,x,a
	return
	end
	
