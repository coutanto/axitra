!*****************************************************************************
!                           AXITRA Force Version
!*                          DIMENSION1
!*
!*         Declaration de variables generales :
!*                constantes, parametres du modele
!*
!*****************************************************************************

module dimension1
   use parameter
   implicit none

   complex(kind=fd)  ::  omega, omega2, a1, xlnf
   real(kind=fd)      ::  kr, kr2, uconv
   logical           ::  ttconv

   integer, allocatable      :: irs(:, :)
   integer, allocatable      :: irc(:), nzr(:), nzrr(:, :), izrr(:, :, :)
   integer, allocatable      :: isc(:), nzs(:), nzss(:, :), izss(:, :, :)
   real(kind=fd), allocatable :: rr(:), hc(:)
   real(kind=fd), allocatable :: cosr(:, :), sinr(:, :)
   real(kind=fd), allocatable :: xs(:), ys(:), zs(:)
   real(kind=fd), allocatable :: xr(:), yr(:), zr(:)
   real(kind=fd), allocatable :: vp(:), vs(:), vp2(:), vs2(:)
   real(kind=fd), allocatable :: rho(:), qp(:), qs(:)

!$OMP THREADPRIVATE(xlnf,omega,omega2,a1,kr,kr2,ttconv)

end module
