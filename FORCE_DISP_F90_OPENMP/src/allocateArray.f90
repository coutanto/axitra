!******************************************************************************
!*                     AXITRA Force Version
!*
!*                         SUBROUTINE ALLOCATEARRAY
!
!
!
module allocateArraym
contains
!******************************************************************************
!                         SUBROUTINE ALLOCATEARRAY1
!
! Allocate array shared between threads
!
! input:
!                nc: number of layer
!                nr: number of receiver
!                ns: number of source
!
!******************************************************************************
   subroutine allocateArray1(nc, nr, ns)

      use dimension1
      use dimension2
      use parameter
      implicit none

      integer :: nc, nr, ns

      allocate (vp(nc), vs(nc), vp2(nc), vs2(nc), rho(nc), qp(nc), qs(nc), hc(nc))
      allocate (xs(ns), ys(ns), zs(ns))
      allocate (xr(nr), yr(nr), zr(nr), cosr(nr, ns), sinr(nr, ns), rr(nr*ns))

      allocate (izrr(nr, nr, nc), irc(nc), nzr(nc), nzrr(nr, nc))
      allocate (izss(ns, ns, nc), isc(nc), nzs(nc), nzss(ns, nc))

      allocate (irs(nr, ns))

      allocate (cff(ns))

   end subroutine
!******************************************************************************
!                         SUBROUTINE ALLOCATEARRAY2
!
! Allocate array private between threads
!
! input:
!                nc: number of layer
!                nr: number of receiver
!                ns: number of source
!
!******************************************************************************
   subroutine allocateArray2(nc, nr, ns, nrs)

      use dimension1
      use dimension2
      use parameter
      implicit none

      integer :: nc, nr, ns, nrs

      allocate (cka(nc))
      allocate (ckb(nc))
      allocate (cka2(nc))
      allocate (ckb2(nc))
      allocate (cnu(nc))
      allocate (cgam(nc))
      allocate (c2(nc))

      allocate (j0(nrs))
      allocate (j1(nrs))
      allocate (k0(nrs))
      allocate (k1(nrs))
      allocate (k2(nrs))
      allocate (k3(nrs))
      allocate (k4(nrs))
      allocate (k5(nrs))

      allocate (rd(nc, 2, 2))
      allocate (ru(nc, 2, 2))
      allocate (td(nc, 2, 2))
      allocate (tu(nc, 2, 2))
      allocate (rdsh(nc))
      allocate (rush(nc))
      allocate (tdsh(nc))
      allocate (tush(nc))
      allocate (me1(nc))
      allocate (me2(nc))

      allocate (nt(nc, 2, 2))
      allocate (mt(nc, 2, 2))
      allocate (ntsh(nc))
      allocate (mtsh(nc))

      allocate (fdo(nc, 2, 2))
      allocate (fup(nc, 2, 2))
      allocate (fupsh(nc))
      allocate (fdosh(nc))

      allocate (su1(ns, 2))
      allocate (sd1(ns, 2))
      allocate (su2(ns, 2))
      allocate (sd2(ns, 2))
      allocate (su3(ns, 2))
      allocate (sd3(ns, 2))
      allocate (su4(ns, 2))
      allocate (sd4(ns, 2))
      allocate (su1sh(ns))
      allocate (sd1sh(ns))
      allocate (su2sh(ns))
      allocate (sd2sh(ns))

      allocate (u(nr, ns, 5))
   end subroutine
end module
