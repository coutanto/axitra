!*****************************************************************************
!                     AXITRA Moment Version
!*
!*                          DIMENSION2
!*
!* 	Declaration de variables utilisees dans le calcul par boucle en kr
!*	externe	 (routines reflect0 a reflect6)
!*
!*****************************************************************************

module dimension2

  implicit none

  complex(kind=8), allocatable :: cka(:),ckb(:),cka2(:),ckb2(:),cnu(:),cgam(:),c2(:),cff(:)
  complex(kind=8), allocatable :: tdsh(:),tush(:),me1(:),me2(:),rdsh(:),rush(:)
  real(kind=8),    allocatable :: fj1(:),k0(:),k1(:),k2(:),k3(:),k4(:),k5(:)
  complex(kind=8), allocatable :: rd(:,:,:),ru(:,:,:),td(:,:,:),tu(:,:,:)
  complex(kind=8), allocatable :: nt(:,:,:),mt(:,:,:),ntsh(:),mtsh(:)
  complex(kind=8), allocatable :: fdo(:,:,:),fup(:,:,:),fupsh(:),fdosh(:)
  complex(kind=8), allocatable :: su1(:,:),sd1(:,:),su2(:,:),sd2(:,:),su3(:,:),sd3(:,:),su4(:,:),sd4(:,:)
  complex(kind=8), allocatable :: su1sh(:),sd1sh(:),su2sh(:),sd2sh(:)
  complex(kind=8), allocatable :: u(:,:,:)


!$OMP THREADPRIVATE(cka,ckb,cka2,ckb2,cnu,cgam,c2)
!$OMP THREADPRIVATE(fj1,k0,k1,k2,k3,k4,k5)
!$OMP THREADPRIVATE(rd,ru,td,tu,rdsh,rush,tdsh,tush,me1,me2)
!$OMP THREADPRIVATE(nt,mt,ntsh,mtsh)
!$OMP THREADPRIVATE(fdo,fup,fupsh,fdosh)
!$OMP THREADPRIVATE(su1,sd1,su2,sd2,su3,sd3,su4,sd4,su1sh,sd1sh,su2sh,sd2sh)
!$OMP THREADPRIVATE(u)

end module
