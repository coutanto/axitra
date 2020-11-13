!*****************************************************************************
!                           AXITRA Force Version
!*                          DIMENSION2
!*
!*         Declaration de variables utilisees dans le calcul par boucle en kr
!*        externe         (routines reflect0 a reflect6)
!*
!*****************************************************************************
module dimension2
   use parameter
   implicit none

   real(kind=fd), allocatable     :: j0(:), j1(:), k0(:), k1(:), k2(:), k3(:), k4(:), k5(:)
   complex(kind=fd), allocatable :: cka(:), ckb(:), cka2(:), ckb2(:), cnu(:), cgam(:), c2(:), cff(:)
   complex(kind=fd), allocatable :: rd(:, :, :), ru(:, :, :), td(:, :, :), tu(:, :, :)
   complex(kind=fd), allocatable :: tdsh(:), tush(:), me1(:), me2(:), rdsh(:), rush(:)
   complex(kind=fd), allocatable :: nt(:, :, :), mt(:, :, :), ntsh(:), mtsh(:)
   complex(kind=fd), allocatable :: fdo(:, :, :), fup(:, :, :), fupsh(:), fdosh(:)
   complex(kind=fd), allocatable :: su1(:, :), sd1(:, :), su2(:, :), sd2(:, :), su3(:, :), sd3(:, :), su4(:, :), sd4(:, :)
   complex(kind=fd), allocatable :: su1sh(:), sd1sh(:), su2sh(:), sd2sh(:)
   complex(kind=fd), allocatable :: u(:, :, :)


!$OMP THREADPRIVATE(cka,ckb,cka2,ckb2,cnu,cgam,c2)
!$OMP THREADPRIVATE(j0,j1,k0,k1,k2,k3,k4,k5)
!$OMP THREADPRIVATE(rd,ru,td,tu,rdsh,rush,tdsh,tush,me1,me2)
!$OMP THREADPRIVATE(nt,mt,ntsh,mtsh)
!$OMP THREADPRIVATE(fdo,fup,fupsh,fdosh)
!$OMP THREADPRIVATE(su1,sd1,su2,sd2,su3,sd3,su4,sd4,su1sh,sd1sh,su2sh,sd2sh)
!$OMP THREADPRIVATE(u)
end module
