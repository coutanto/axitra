!
!******************************************************************************
!
!  read all axitra input files used for computing Green's functions
!   and format them for an intregation into python
!   All data needed are serialized as follow
!
!
!******************************************************************************

program axitra2py

   use dimension1
   use dimension2
   use parameter
   use initdatam
   use allocatearraym


   implicit none
! Local
   character(len=20)    :: sourcefile, statfile, header,arg
   integer              :: ic, ir, is, nfreq, ikmax, ncp, iklast, jf, ik, lastik
   integer              :: nrs ! number of receiver radial distance
   integer              :: ncr ! number of layer containing a receiver
   integer              :: ncs ! number of layer containing a source

   integer              :: nr, ns, nc, narg, out3
   real(kind=fd)         :: dfreq, freq, pil
   logical              :: latlon, freesurface
   logical, allocatable :: tconv(:, :)
   real(kind=fd)         :: rw, aw, phi, zom, tl, xl, rmax, vlim, vmean
   namelist/input/nfreq, tl, aw, xl, ikmax, latlon, freesurface, sourcefile, statfile

   integer, allocatable :: index(:),rindex(:)

#include "version.h"
   write(0,*) 'running axitra2py '//VERSION

!
! read header if any
!
   narg=iargc()
   if (narg>0) then
     call getarg(1,arg)
     write(header,"('axi_',A)") trim(arg)
   else
     header='axi'
   endif
!++++++++++
!           Read input parameter from <header>.data
!           <header>.head is used later to know the exact number of frequency
!           actually computed
!
!++++++++++
   open (in1, form='formatted', file=trim(header)//'.data')
   out3=6
   open (out3,  form='formatted', file='axi.datapy')

   if (narg>0) then
      write(out3,*) trim(arg)
   else
      write(out3,*) 'no_suffix'
   endif

! count number of layer
   read (in1, input, err=100)
   write (out3, *) nfreq, tl, aw, xl, ikmax, latlon, freesurface
   write (out3, *) sourcefile
   write (out3, *) statfile
   nc=0
   do while(.true.)
     read(in1,*,end=91)
     nc=nc+1
   end do
91 rewind(in1)

   open (in2, form='formatted', file=sourcefile)
   open (in3, form='formatted', file=statfile)

! count number of sources
   ns=0
   do while(.true.)
      read(in2,*,end=92)
      ns=ns+1
   end do
92 rewind(in2)

! count number of receiver
   nr=0
   do while(.true.)
     read(in3,*,end=93)
     nr=nr+1
   end do
93 rewind(in3)
   call allocateArray1(nc, nr, ns)

!
! read velocity model
!
   read (in1, input)
   write(out3, *) nc
   do ic = 1, nc
      read (in1, *) hc(ic), vp(ic), vs(ic), rho(ic), qp(ic), qs(ic)
      write (out3, *) hc(ic), vp(ic), vs(ic), rho(ic), qp(ic), qs(ic)
   enddo
!++++++++++++
!        read stations and source coordinates
!++++++++++++

   write(out3, *) ns
   allocate(index(ns))
   do is = 1, ns
      read (in2, *) index(is), xs(is), ys(is), zs(is)
      write(out3, *)   index(is), xs(is), ys(is), zs(is)
   enddo
   write(out3,*) nr
   allocate(rindex(nr))
   do ir = 1, nr
      read (in3, *) rindex(ir), xr(ir), yr(ir), zr(ir)
      write(out3, *)   rindex(ir), xr(ir), yr(ir), zr(ir)
   enddo

!++++++++++
!           INITIALISATIONS
!++++++++++
   rewind(in2)
   rewind(in3)
   call initdata(latlon, nr, ns, nc, ncr, ncs,nrs,rmax)

! compute xl and or tl if not supplied
   if (xl<=0.d0 .or. tl<=0.d0) call estimateParameter(xl,tl,rmax,vp,vs,nc)
   write(out3, *) xl,tl

   close(out3)
   close(in1)
   close(in2)

   call exit(0)
100 call exit(1)
end

!
! estimate tl and/or xl parameter in case they were not given by user
!

subroutine estimateParameter(xl,tl,rmax,vp,vs,nc)
use parameter
implicit none
real(kind=fd) :: xl,tl,rmax,vp(nc),vs(nc)
integer      :: nc

integer      :: i
real(kind=fd) :: vmax,vmin

vmax=vp(1)
vmin=vs(1)
do i=2,nc
   vmax=max(vmax,vp(i))
   vmin=min(vmin,vs(i))
enddo

! estimate a duration
! we need at least the time needed to travel
! along the largest distance at the lowest velocity
vmin=0.8d0*vmin
if (tl<=0.d0) then
    tl=rmax/vmin
    tl=(int(tl/5.d0)+1)*5.d0
write(6,*) 'duration tl set automatically to ',tl,'sec'
endif

! estimate a radial periodicity
! source are located with a xl periodicity
! To avoid the effect of periodic sources during the
! tl duration, we need:
if (xl<=0.d0) then
    xl = 1.2d0*rmax + vmax*tl
    xl=(int(xl/10.d0)+1)*10.d0
    write(6,*) 'periodicity xl set automatically to ',xl,'(k)m'
endif
end subroutine
