!**********************************************************************
!*                  PARAMETER
!*
! *    ikmin: minimum iteration number for summation
! *    fref: reference frequency for attenuation
!**********************************************************************

module parameter
   use iso_fortran_env
   implicit none

! size for kind declaration, fs= single precision, fd=double precision
! uses iso_fortran_env, can be replace by any sizeof if your compiler is not Fortran2003 compliant
   integer :: fs,fd
   parameter(fd=real64, fs=real32)

   complex(kind=fd)  ::  ai
   parameter(ai=(0., 1.))
   real(kind=fd)      ::  pi, pi2
   parameter(pi=3.14159265359d0, pi2=6.28318530718d0)
   integer :: in1, in2, in3, out, out2
   parameter(in1=10, in2=11, in3=12, out=13, out2=14)

   integer :: ikmin
   parameter(ikmin=100)

! explim= exponent lower limit
! needs to be adjuest when running with single precision float number
   real(kind=fd) :: explim, elim
   parameter(explim=-300., elim=1.d-300)

! bessel functions are stored up to nkmax wavenumber iterations.
! the larger, the faster, but needs more memory
   integer :: nkmax
   parameter(nkmax=2000) ! not used yet

! convergence relative error
   real(kind=fd) :: rerr
   parameter(rerr=1.e-4)

! reference frequency for attenuation
   real(kind=fd) :: fref
   parameter(fref=1.)

contains

function doubleEquality(x1,x2) result(t)
use ISO_FORTRAN_ENV
implicit none

real(kind=real64) :: x1,x2
logical :: t

t=(abs(x1-x2)<=epsilon(x1))
end function
end module
