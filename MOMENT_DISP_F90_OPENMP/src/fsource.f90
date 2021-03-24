!**********************************************************
!	FSOURCE
!
!	Different kind of source function defined in the
!	frequency domain
!	Source function for dislocation (step, ramp, haskell)
!	are normalized so that in far-field the low-frequency
!	level is proportionnal to the seismic moment with a factor
!	equals to: Rad/(4.PI.rho.beta^3) * 1./r
!	where Rad= Radiation coefficient with (possible)
!	free surface effect
!
!	input:
!		type 	-> 	see below
!		omega	->	angular frequency
!		t0,t1	->	time constant when needed
!		dt	->	time step
!               
!
!  ATTENTION
!        Spectra are defined with respect to the FFT normalization 
!        used in the code (1/N for inverse Fourier Transform, N= number of points)
!        E.G. For a dirac source time function, the spectrum is set to 
!        to modulus = 1.
!        The spectrum is later multiplied by 1/nt
!        
!**********************************************************
module fsourcem

contains
function	fsource (type, t0, omega, t1, dt)
    use parameter
    implicit 	none
    integer          ::	type
    real(kind=fd)    :: dt,t0,t1
    real(kind=fd)    ::  rand
    real(kind=fd)    :: uur,uui,trise,trupt
    complex(kind=fd) :: fsource,uu,uex,uxx,omega,shx
    real(kind=fd),parameter    :: two=2.d0

! TYPE=0               Source time function = unit dirac
!                      For a shear dislocation
!                      source slip is a dirac,
!                      i.e slip and return back to initial position
!                      
    if (type.eq.0) then
           fsource=1.d0
    endif
 
! TYPE=1        Source time function = unit  Ricker
!               (Double derivative of a Gaussian pulse)
!
    if (type.eq.1) then
          uu=t0*omega/two
          uu=exp(-uu*uu)
          uu=uu*t0*sqrt(pi)
          uu=uu*omega*omega*t0*t0/two
          fsource=uu
    endif
 
! TYPE=2        Source time function =  smooth non causal unit step
!		(ArcTang function with symetrical point at t=0)

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
    endif

! TYPE=7        Source = Heaviside / unit step
    if (type.eq.7) then
	  uu=1./ai/omega
          fsource= uu
    endif

! TYPE=4        Source = integral of a triangle 
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
 
! TYPE=5        Source = causal ramp
!               rise time T=t0
    if (type.eq.5) then
	  trise=t0
          uu=ai*omega*trise
          uu=(1.-exp(-uu))/uu
          fsource=uu/(ai*omega)
    endif
	
! TYPE=6,8        Source = kind of Haskell, trapezoide
!	1 st time cste, rise time: trise
!	2 nd time cste, rupture duration: trupt
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


    return
end function
end module
