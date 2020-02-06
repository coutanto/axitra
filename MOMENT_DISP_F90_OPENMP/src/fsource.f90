!**********************************************************
!	FSOURCE
!
!	Different kind of source function defined in the
!	frequency domain
!	Source function for dislocation (step, ramp, haskell)
!	are normalized so that in far-field the low-frequency
!	level is proportionnal to the seismic moment with a factor
!	equals to: Rad/(4.PI.rho.beta^3) * 1./r
!	where Rad= Radiation coefficient with (possibly)
!	free surface effect
!
!	input:
!		type 	-> 	see below
!		omega	->	angular frequency
!		t0,t1	->	time constant when needed
!		dt	->	sampling rate
!**********************************************************
module fsourcem

contains
function	fsource (type, t0, omega, t1, dt)
use parameter
	implicit 	none
	integer		type
	real(kind=8)		dt,t0,t1
	real rand
	real(kind=8)    :: uur,uui,trise,trupt
	complex(kind=8) :: fsource,uu,uex,uxx,omega,shx

!       write(0,*) type, omega, t0, t1, dt
! TYPE=0               Source = Dirac en deplacement
        if (type.eq.0) then
          fsource=1
        endif
 
! TYPE=1        Source = Ricker en deplacement
        if (type.eq.1) then
          uu=omega*t0
          uu=uu*uu/pi2/pi2
          uu=exp(-uu)
          uu=omega*omega*uu*dt
	  fsource= uu
        endif
 
! TYPE=2        Source = step en deplacement
!		2 steps possibles (1) real=1/(ai*omega)
!				  (2) bouchon s
        if (type.eq.2.or.type.eq.9) then
!         if (zabs(omega*pi*t0/2.).lt.700) then
          if (zabs(omega*pi*t0/2.).lt.20) then
          	shx=zexp(omega*pi*t0/2.)	!Bouchon s
                shx=1./(shx-1./shx)
          else
                shx=0.
          endif
          uu=-ai*t0*pi*shx
          fsource= uu
!         write(0,*) omega,t0,pi,omega*pi*t0/2.,shx
        endif
! TYPE=7        Source = step en deplacement
	if (type.eq.7) then
	  uu=1./ai/omega
          fsource= uu
	endif

! TYPE=3        Source = fichier  axi.sou
!               sismogramme dans l unite choisie dans le fichier
        if (type.eq.3) then
          read(16,*) uur,uui
          fsource=cmplx(uur,uui)
        endif
 
! TYPE=4        Source = triangle en deplacement
        if (type.eq.4) then
	  trise=t0
	  trupt=t0
	  uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu		! ramp
          uxx=ai*omega*trupt/2.		! finite fault
          uex=exp(uxx)
          uxx=(uex-1./uex)/uxx/2.
          fsource=uu*uxx/(ai*omega)
        endif
 
! TYPE=5        Source = rampe causale
!               rise time T=t0
        if (type.eq.5) then
	  trise=t0
          uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu
          fsource=uu/(ai*omega)
        endif
	
! TYPE=6,8        Source = modele d haskell, trapezoide
!	1 ere cste de temps rise time: riset
!	2 eme cste de temps, duree de la rupture
!	  trupt = Length/2/rupt_velocity (Haskell)
        if ((type.eq.6).or.(type.eq.8)) then
	  trise=t0
	  trupt=t1
	  uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu		! ramp
          uxx=ai*omega*trupt/2.		! finite fault
          uex=exp(uxx)
          uxx=(uex-1./uex)/uxx/2.
          fsource=uu*uxx/(ai*omega)
        endif
! TYPE 9 bruit blanc
	if (type.eq.9) then
          fsource=fsource*exp(ai*rand(0)*2*pi)
	endif
    
        return
	end
end module
