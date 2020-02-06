!*****************************************************************************
!                           AXITRA Force Version
!*                          DIMENSION1
!*
!*         Declaration de variables generales :
!*                constantes, parametres du modele
!*
!*****************************************************************************

module dimension1

   implicit none

   complex(kind=8)  ::  omega, omega2, a1, xlnf
   real(kind=8)      ::  kr, kr2, uconv
   logical           ::  ttconv

   integer, allocatable      :: irs(:, :)
   integer, allocatable      :: irc(:), nzr(:), nzrr(:, :), izrr(:, :, :)
   integer, allocatable      :: isc(:), nzs(:), nzss(:, :), izss(:, :, :)
   real(kind=8), allocatable :: rr(:), hc(:)
   real(kind=8), allocatable :: cosr(:, :), sinr(:, :)
   real(kind=8), allocatable :: xs(:), ys(:), zs(:)
   real(kind=8), allocatable :: xr(:), yr(:), zr(:)
   real(kind=8), allocatable :: vp(:), vs(:), vp2(:), vs2(:)
   real(kind=8), allocatable :: rho(:), qp(:), qs(:)

!$OMP THREADPRIVATE(xlnf,omega,omega2,a1,kr,kr2,ttconv)

end module
