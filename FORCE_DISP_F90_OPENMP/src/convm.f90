!******************************************************************************
!  PROGRAMME CONVM
!
!  FORCE version
!
! Convolution and output using output() routine
! All computations done in double float format(64bits)
! but output is converted to 32bits
!
! nsp = max number of sources
! ncp = max number of layers
! nrp = max number of receiver
! ntp = max length(in points) for seismogram
!******************************************************************************
#define NSTYPE 3
program convm

   use parameter
   use fsourcem
   implicit none

   integer ntp, nrtp, icc, ics, ic, base

   character filein*10, sourcefile*20, statfile*20, chan(3)*2

   integer:: jf, ir, is, it, nc, ns, nr, nfreq, ikmax, mm, nt, io, index, indexin
   real(kind=8)    ::  tl, xl, uconv, hh, zsc, dfreq, freq, aw, ck, xmm, xref, yref,  &
                       lat, long, t0, t1, hanning, xphi, dt0, pas
   complex(kind=8) ::  omega, uxf(NSTYPE), uyf(NSTYPE), uzf(NSTYPE), deriv, us, uux, uuy, uuz, cc, fs
   logical         ::  latlon,freesurface
   integer, allocatable         :: iwk(:), isc(:), rindex(:)
   real(kind=8), allocatable    :: hc(:), vp(:), vs(:), rho(:), delay(:), xr(:), yr(:), &
                                   zr(:), a(:, :), qp(:), qs(:), xs(:), ys(:), zs(:), amp(:)
   real(kind=4), allocatable    :: sx(:), sy(:), sz(:)
   real(kind=4)                 :: spas
   complex(kind=8), allocatable :: ux(:, :), uy(:, :), uz(:, :)

   namelist/input/nc, nfreq, tl, aw, nr, ns, xl, ikmax, latlon, freesurface, sourcefile, statfile

   dt0 = 0.
   chan(1) = 'X'
   chan(2) = 'Y'
   chan(3) = 'Z'

! read binary results from axitra
! We assume here that record length is given in byte.
! For intel compiler, it means using "assume byterecl" option
! record length is 3 x (complex double precision) = 3 x 2 x 8 bytes
   open (12, access='direct', recl=NSTYPE*2*8, form='unformatted', file='axi.res')

! read input file to axitra
   open (10, form='formatted', file='axi.data')
   read (10, input)

! allocate space knowing number of sources, receivers and layers
   allocate (hc(nc), vp(nc), vs(nc), qp(nc), qs(nc), rho(nc))
   allocate (delay(ns), xs(ns), ys(ns), zs(ns), a(NSTYPE , ns), amp(ns),isc(ns))
   allocate (xr(nr), yr(nr), zr(nr), rindex(nr))

! other input files
   open (13, form='formatted',file=sourcefile)
   open (14, form='formatted', file=statfile)
   open (15, form='formatted', file='axi.hist')


! + +++++++++++
! read source function type from keyboard
! + +++++++++++

   write (6, *) 'source function?'
   write (6, *) '0 : Dirac'
   write (6, *) '1 : Ricker'
   write (6, *) '2 : step '
   write (6, *) '3 : function stored in file <axi.sou>'
   write (6, *) '4 : triangle'
   write (6, *) '5 : ramp'
   write (6, *) '6 : not used....'
   write (6, *) '7 : True step (watch high frequencies cutoff!!)'
   read (5, *) ics

   if (ics .eq. 1) then
      write (6, *) "Ricker pseudo period?"
      read (5, *) t0
   endif
   if (ics .eq. 2) then
      write (6, *) "source duration?"
      read (5, *) t0
   endif
   if (ics .eq. 3) then
      open (16, file='axi.sou',access='stream')
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
      read (13, *) index, xs(is), ys(is), zs(is)
      indexin = -1
      rewind (15)
      do while (indexin .ne. index)
!  read: index, fx , fy, fz, amplitude, time_delay
         read (15, *) indexin, a(1,is), a(2,is), a(3,is), amp(is), delay(is)
      enddo
      delay(is) = delay(is) + dt0
   enddo

   do ir = 1, nr
      read (14, *) rindex(ir), xr(ir), yr(ir), zr(ir)
! reference: x = north; y = east
   enddo
!
   if (latlon) then
      call ll2km(xr, yr, nr, xs, ys, ns)
   endif

!++++++++++++
!        Parametres sources
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

      a(1,is)=a(1,is)*amp(is)
      a(2,is)=a(2,is)*amp(is)
      a(3,is)=a(3,is)*amp(is)
   enddo

!++++++++++++
!        Lecture fonctions de transfert
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

   dfreq = 1./tl
   pas = tl/nt
   spas=sngl(pas)
   aw = -pi*aw/tl

! loop over frequencies
   open (10, form='formatted', file='axi.head')
   do jf = 1, nfreq
      freq = float(jf - 1)/tl
      omega = cmplx(pi2*freq, aw)
      deriv = (ai*omega)**icc

      read (10, *, end=2000)
	  fs=fsource(ics, t0, omega, t1, pas)
      do is = 1, ns
         xphi = atan2(dble(yr(1) - ys(is)), dble(xr(1) - xs(is)))
!         rvel(is) = rpvel*vs(isc(is))/(1.-rpvel*cos(xphi - strike(is)))
!         if (ics .ne. 8) then
!            t1 = length(is)/rvel(is)
!         endif

         us = fs*deriv*exp(-ai*omega*delay(is))

         do ir = 1, nr
            base = (ir - 1)*3 + (is - 1)*3*nr + (jf - 1)*3*nr*ns
            read (12, rec=base + 1, err=2000) (uxf(it), it=1, NSTYPE)
            read (12, rec=base + 2, err=2000) (uyf(it), it=1, NSTYPE)
            read (12, rec=base + 3, err=2000) (uzf(it), it=1, NSTYPE)
            uux = 0.
            uuy = 0.
            uuz = 0.
            do it = 1, NSTYPE
               uux = uux + uxf(it)*a(it, is)
               uuy = uuy + uyf(it)*a(it, is)
               uuz = uuz + uzf(it)*a(it, is)
               if (isnan(real(uux)) .or. isnan(real(uuy)) .or. isnan(real(uuz))) write(0,*) 'real part uux,uuy or uuz is NaN'
               if (isnan(imag(uux)) .or. isnan(imag(uuy)) .or. isnan(imag(uuz))) write(0,*) 'imag part uux,uuy or uuz is NaN'
            enddo

            ux(jf, ir) = ux(jf, ir) + uux*us
            uy(jf, ir) = uy(jf, ir) + uuy*us
            uz(jf, ir) = uz(jf, ir) + uuz*us

         enddo !fin boucle recepteur
      enddo !fin boucle source
   enddo
2000 continue
   if (jf < nfreq) then
      nfreq = jf - 1
   endif

!++++++++++++
!                Calcul des sismogrammes
!++++++++++++

   do ir = 1, nr

!                on complete le spectre pour les hautes frequences
!               avec inversion du signe de la partie imaginaire pour
!               la FFT inverse qui n existe pas avec fft2cd
      do jf = nt + 2 - nfreq, nt
         ux(jf, ir) = conjg(ux(nt + 2 - jf, ir))
         uy(jf, ir) = conjg(uy(nt + 2 - jf, ir))
         uz(jf, ir) = conjg(uz(nt + 2 - jf, ir))
      enddo

      call fft2cd(ux(1, ir), mm, iwk)
      call fft2cd(uy(1, ir), mm, iwk)
      call fft2cd(uz(1, ir), mm, iwk)

      do it = 1, nt
         ck = float(it - 1)/nt
         cc = exp(-aw*tl*ck)/tl
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


real function hanning(i, max, perc)
   use parameter
   implicit none

   integer i, max
   real(kind=8) xi, xm, perc
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
