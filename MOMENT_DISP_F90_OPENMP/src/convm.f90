
!   PROGRAMME CONVM
!
!   Moment version
!
! Convolution and output using output() routine
! All computations done in double float format(64bits)
! but output is converted to 32bits
!
! convm read input files: <header>.data, <header>.res, <header>.hist, <header>.sou
! where <header>='axi' or 'axi_???' and ??? is given as first argument to convm
!
! nsp = max number of sources
! ncp = max number of layers
! nrp = max number of receiver
! ntp = max length(in points) for seismogram
!******************************************************************************

! number of elementary source, 6 for moment, 3 for forces
#define NSTYPE 6 
program convm

   use parameter
   use fsourcem
   implicit none

   integer ntp, nrtp, icc, ics, ic, base, narg

   character header*10, sourcefile*20, statfile*20, chan(3)*2,arg*10

   integer:: jf, ir, is, it, nc, ns, nr, nfreq, ikmax, mm, nt, io, indexin
   real(kind=fd)    ::  tl, xl, uconv, hh, zsc, dfreq, freq, aw, ck, xmm, xref, yref, &
                        lat, long, t0, t1, hanning, pas, xphi, dt0, rfsou, stime
   complex(kind=fd) ::  omega, uxf(NSTYPE), uyf(NSTYPE), uzf(NSTYPE), deriv, us, uux, uuy, uuz, cc, freqs
   logical                   :: latlon,freesurface
   integer, allocatable      :: iwk(:), isc(:), rindex(:), sindex(:)
   real(kind=fd), allocatable :: hc(:), vp(:), vs(:), rho(:), delay(:), xr(:), yr(:), zr(:), a(:, :), qp(:), qs(:), &
                                mu(:), strike(:), dip(:), rake(:), disp(:), xs(:), ys(:), zs(:), width(:), length(:)
   real(kind=fs), allocatable :: sx(:), sy(:), sz(:)
   real(kind=fs)              :: spas
   complex(kind=fd), allocatable :: ux(:, :), uy(:, :), uz(:, :), fsou(:)

   namelist/input/ nfreq, tl, aw, xl, ikmax, latlon, freesurface, sourcefile, statfile

   dt0 = 0.
   chan(1) = 'X'
   chan(2) = 'Y'
   chan(3) = 'Z'

! header if any
   narg=iargc()
   if (narg>0) then
     call getarg(1,arg)
     write(header,"('axi_',A)") trim(arg)
   else
     header='axi'
   endif

! read binary results from axitra
! We assume here that record length is given in byte.
! For intel compiler, it means using "assume byterecl" option
! record length is 6 x (complex double precision) = 6 x 2 x 8 bytes
   open (12, access='direct', recl=NSTYPE*2*8, form='unformatted', file=trim(header)//'.res')

! read input file to axitra
   open (10, form='formatted', file=trim(header)//'.data')
! count number of layer
   read (10, input)
   nc=0
   do while(.true.)
     read(10,*,end=91)
     nc=nc+1
   end do
91 rewind(10)
   read (10, input)

   open (13, form='formatted', file=sourcefile)
   open (14, form='formatted', file=statfile)
! count number of sources
   ns=0
   do while(.true.)
     read(13,*,end=92)
     ns=ns+1
   end do
92 rewind(13)
! count number of receiver
   nr=0
   do while(.true.)
     read(14,*,end=93)
     nr=nr+1
   end do
93 rewind(14)

! allocate space knowing number of sources, receivers and layers
   allocate (hc(nc), vp(nc), vs(nc), qp(nc), qs(nc), rho(nc))
   allocate (delay(ns), mu(ns), strike(ns), rake(ns), dip(ns), disp(ns))
   allocate (xs(ns), ys(ns), zs(ns), isc(ns), a(NSTYPE, ns), width(ns), length(ns))
   allocate (xr(nr), yr(nr), zr(nr), rindex(nr),sindex(ns))

! other input files

   open (15, form='formatted', file=trim(header)//'.hist')

! + +++++++++++
! read source function type from keyboard
! + +++++++++++

   write (6, *) 'Select the Source time function'
   write (6,*)  'Add +10 to compute only the time function (no convolution)'
   write (6, *) '0 : Dirac'
   write (6, *) '1 : Ricker'
   write (6, *) '2 : step '
   write (6, *) '3 : function stored in file <axi.sou>'
   write (6, *) '4 : triangle'
   write (6, *) '5 : ramp'
   write (6, *) '6 : not used....'
   write (6, *) '7 : True step (watch high frequencies cutoff!!)'
   write (6, *) '8 : Trapezoid'
   read (5, *) ics
   
   ! test wether we want only the source time function
   if (ics>=10) then
       ics=ics-10
       stime = 0.d0
   else
       stime = 1.d0
   endif

   if (ics .eq. 8) then
      write (6, *) "rise time t0 ?"
      read (5, *) t0
      write (6, *) "source process time t1 ? (total time = t1+t0)"
      read (5, *) t1
      t1 = t1 + t0
   endif

   if (ics .eq. 1) then
      write (6, *) "Ricker pseudo period?"
      read (5, *) t0
   endif
   if (ics .eq. 2) then
      write (6, *) "source duration?"
      read (5, *) t0
   endif

   if (ics .eq. 4) then
      write (6, *) "triangle width ?"
      read (5, *) t0
   endif
   if (ics .eq. 5) then
      write (6, *) "rise time?"
      read (5, *) t0
   endif

   write (6, *) 'output is ?'
   write (6, *) '1: displacement (m)'
   write (6, *) '2: velocity (m/s)'
   write (6, *) '3: acceleration (m/s/s)'
   read (5, *) icc
   icc = icc - 1

!++++++++++++
!                 medium, source and station
!++++++++++++
   do ic = 1, nc
      read (10, *) hc(ic), vp(ic), vs(ic), rho(ic), qp(ic), qs(ic)
   enddo
   close (10)

   do is = 1, ns
      read (13, *) sindex(is), xs(is), ys(is), zs(is)
   enddo
! sort by increasing depth because Green's functions are sorted by source increasing depth
   call sortByDepth(sindex,xs,ys,zs,ns)

   do is = 1,ns
      indexin = -1
      rewind (15)
      do while (indexin .ne. sindex(is))
         read (15, *) indexin, disp(is), strike(is), dip(is), rake(is), width(is), length(is), delay(is)
      enddo
      delay(is) = delay(is) + dt0
   enddo

   do ir = 1, nr
      read (14, *) rindex(ir), xr(ir), yr(ir), zr(ir)
! reference: x = north; y = east
   enddo
! sort by increasing depth because Green's functions are sorted by receiver increasing depth
   call sortByDepth(rindex,xr,yr,zr,nr)
!
   if (latlon) then
      call ll2km(xr, yr, nr, xs, ys, ns)
   endif

!++++++++++++
!        source Parameters
!++++++++++++
!     conversion depth -> thickness
   if (hc(1) .eq. 0.) then
      do ic = 1, nc - 1
         hc(ic) = hc(ic + 1) - hc(ic)
      enddo
   endif

   do is = 1, ns
      hh = 0.
      isc(is) = 1
      zsc = zs(is)
      do ic = 1, nc - 1
         hh = hc(ic)
         if (zsc .gt. hh) then
            zsc = zsc - hh
            isc(is) = ic + 1
         else
            exit
         endif
      enddo

!     rvel(is)=rpvel*vs(isc(is))*5.         !5 a 1 suivant direction
      mu(is) = vs(isc(is))*vs(isc(is))*rho(isc(is))
      call cmoment(mu(is), strike(is), dip(is), rake(is), disp(is), width(is)*length(is), a(1, is))
   enddo

!++++++++++++
! read spatial periodicity and time duration
!++++++++++++
   open (10, form='formatted', file=trim(header)//'.head')
   read (10,*) xl,tl

!++++++++++++
! Initialize time and dimension parameters
!++++++++++++
   xmm = log(real(nfreq))/log(2.)
   mm = int(xmm) + 1
   if (xmm-mm+1 >0) mm=mm+1
   nt = 2**mm
   allocate (iwk(nt), sx(nt), sy(nt), sz(nt))
   allocate (ux(nt, nr), uy(nt, nr), uz(nt, nr))
   ux = 0.d0
   uy = 0.d0
   uz = 0.d0
   aw = -pi*aw/tl
   dfreq = 1./tl
   pas = tl/nt
   spas=sngl(pas)

! if source time function is supplied, read it
   if (ics .eq. 3) then
     write(0,*) 'convms try reading ',nt,'time points from file ',trim(header)//'.sou'
     open (16, form='formatted', file=trim(header)//'.sou')
     allocate(fsou(nt))
     fsou=0.d0
     do it=1,nt
       read(16,*, end=1960) rfsou
       ck = float(it - 1)/nt
       cc = exp(aw*tl*ck)
       fsou(it)=rfsou*cc
     enddo
     ! direct FFT normalization
     ! Since the same call is used for inverse transform
     ! needs conjugate
     call fft2cd(fsou, mm, iwk)
     fsou=conjg(fsou)

1960 continue
   endif


! loop over frequencies

   do jf = 1, nfreq
      freq = float(jf - 1)/tl
      omega = cmplx(pi2*freq, aw)
      deriv = (ai*omega)**icc

      read (10, *, end=2000)
      if (ics.eq.3) then
        freqs=fsou(jf)
      else
        freqs=fsource(ics, t0, omega, t1, pas)
	!if (ics==0) freqs = freqs / nt
      endif

      do is = 1, ns ! loop over source
         xphi = atan2(dble(yr(1) - ys(is)), dble(xr(1) - xs(is)))
!         rvel(is) = rpvel*vs(isc(is))/(1.-rpvel*cos(xphi - strike(is)))
!         if (ics .ne. 8) then
!            t1 = length(is)/rvel(is)
!         endif
         us = freqs*deriv*exp(-ai*omega*delay(is))
         do ir = 1, nr ! loop over receiver
            base = (ir - 1)*3 + (is - 1)*3*nr + (jf - 1)*3*nr*ns
            read (12, rec=base + 1, err=2000) (uxf(it), it=1, NSTYPE)
            read (12, rec=base + 2, err=2000) (uyf(it), it=1, NSTYPE)
            read (12, rec=base + 3, err=2000) (uzf(it), it=1, NSTYPE)
            uux = 0.
            uuy = 0.
            uuz = 0.
            do it = 1, NSTYPE
               uux = uux + (uxf(it)*stime + (1.d0-stime))*a(it, is)
               uuy = uuy + (uyf(it)*stime + (1.d0-stime))*a(it, is)
               uuz = uuz + (uzf(it)*stime + (1.d0-stime))*a(it, is)
               if (isnan(real(uux)) .or. isnan(real(uuy)) .or. isnan(real(uuz))) write(0,*) 'real part uux,uuy or uuz is NaN'
               if (isnan(imag(uux)) .or. isnan(imag(uuy)) .or. isnan(imag(uuz))) write(0,*) 'imag part uux,uuy or uuz is NaN'
            enddo
            ux(jf, ir) = ux(jf, ir) + uux*us
            uy(jf, ir) = uy(jf, ir) + uuy*us
            uz(jf, ir) = uz(jf, ir) + uuz*us

         enddo !end receiver loop
      enddo !end source loop
   enddo
2000 continue
   if (jf < nfreq) then
      nfreq = jf - 1
   endif

!++++++++++++
!                compute time domain traces
!++++++++++++

   do ir = 1, nr

!               fill-in symetric part of the spectrum
!               reverse imaginery part to perform inverse FFT
      do jf = nt + 2 - nfreq, nt
         ux(jf, ir) = conjg(ux(nt + 2 - jf, ir))
         uy(jf, ir) = conjg(uy(nt + 2 - jf, ir))
         uz(jf, ir) = conjg(uz(nt + 2 - jf, ir))
      enddo

! inverse fourier transform
      call fft2cd(ux(1, ir), mm, iwk)
      call fft2cd(uy(1, ir), mm, iwk)
      call fft2cd(uz(1, ir), mm, iwk)

! correct for time domain effect of Green's function imaginery part (aw)
! and apply Fourier transform normalization (1/nt)
      do it = 1, nt
         ck = float(it - 1)/nt
         cc = exp(-aw*tl*ck)
         sx(it) = real(ux(it, ir)*cc)
         sy(it) = real(uy(it, ir)*cc)
         sz(it) = real(uz(it, ir)*cc)
      enddo

      call output(sx, rindex(ir), spas, nt, 11, chan(1))
      call output(sy, rindex(ir), spas, nt, 11, chan(2))
      call output(sz, rindex(ir), spas, nt, 11, chan(3))

   enddo

   stop
end

subroutine cmoment(mu, strike, dip, rake, disp, surf, a)
   use parameter
   implicit none

   real(kind=fd) xmoment, mu, strike, dip, rake, disp, surf, a(6), sd, &
      cd, sp, cp, s2p, s2d, c2p, c2d, x1, x2, x3, x4, x5, x6, cl, sl

   if (surf .eq. 0.) then
      xmoment = disp
   else
      xmoment = mu*disp*surf
   endif
!        xmoment=1.e30
   write (6, *) "moment (Nm):", xmoment
   write (6, *) "moment (Dyne.cm):", xmoment*1.e7
   strike = strike*pi/180.
   dip = dip*pi/180.
   rake = rake*pi/180.
   sd = sin(dip)
   cd = cos(dip)
   sp = sin(strike)
   cp = cos(strike)
   sl = sin(rake)
   cl = cos(rake)
   s2p = 2.*sp*cp
   s2d = 2.*sd*cd
   c2p = cp*cp - sp*sp
   c2d = cd*cd - sd*sd

!       Coefficient pour les sources Mxx,Mxy,Mxz,Myy,Myz,Mzz
   x1 = -(sd*cl*s2p+s2d*sl*sp*sp)*xmoment
   x2 = (sd*cl*c2p+s2d*sl*s2p/2.)*xmoment
   x3 = -(cd*cl*cp + c2d*sl*sp)*xmoment
   x4 = (sd*cl*s2p-s2d*sl*cp*cp)*xmoment
   x5 = -(cd*cl*sp - c2d*sl*cp)*xmoment
   x6 = (s2d*sl)*xmoment

!       Coefficient pour les sources bis (5 dislocations elementaires
!       et une source isotrope)
   a(1) = x2
   a(2) = x3
   a(3) = -x5
   a(4) = (-2.*x1 + x4 + x6)/3.
   a(5) = (x1 - 2*x4 + x6)/3.
!       a(6) = (x1+x4+x6)/3.
   a(6) = 0.
!   If surf=-1. it is an explosion
   if (surf .eq. -1.) then
      xmoment = disp
      a(1) = 0.
      a(2) = 0.
      a(3) = 0.
      a(4) = 0.
      a(5) = 0.
      a(6) = xmoment
!   If surf=-2. it is a tensile crack whose normal
!   should be given
   else if (surf .eq. -2.) then
   endif

   return
end

real function hanning(i, max, perc)
   use parameter
   implicit none

   integer i, max
   real(kind=fd) xi, xm, perc
   xi = i
   xm = max
   if ((xm - xi)/xm .le. perc) then
      hanning = (cos((1 - (xm - xi)/xm/perc)*pi) + 1)/2.
!        else if (xi/xm.le.perc) then
!                hanning=(cos((1-(xi)/xm/perc)*pi)+1)/2.
!        else if (xi.le.10) then
!                hanning=0.
   else
      hanning = 1.
   endif
end
subroutine sortByDepth(rindex,xr,yr,zr,nr)
  implicit none
  real(kind=8) :: xr(*),yr(*),zr(*),tmp
  integer      :: nr,ir,jr,itmp,rindex(*)
!++++++++++++
!         sort stations according to increasing depth
!         we do that whatever the ordering was in the station file
!++++++++++++
   do ir = 1, nr - 1
      do jr = ir, nr
         if (zr(ir) .gt. zr(jr)) then
            tmp = xr(ir)
            xr(ir) = xr(jr)
            xr(jr) = tmp
            tmp = yr(ir)
            yr(ir) = yr(jr)
            yr(jr) = tmp
            tmp = zr(ir)
            zr(ir) = zr(jr)
            zr(jr) = tmp
            itmp = rindex(ir)
            rindex(ir) = rindex(jr)
            rindex(jr) = itmp
         endif
      enddo
   enddo
end subroutine
