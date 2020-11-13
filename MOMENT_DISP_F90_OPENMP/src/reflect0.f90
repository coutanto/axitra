!******************************************************************************
!                     AXITRA Moment Version
!*
!*                     SUBROUTINE REFLECT0
!*
!*    Calcul de coefficients dependant du nombre d onde kr et de fonctions
!*    de Bessel.
!*
!******************************************************************************
module reflect0m
use parameter
contains
subroutine initc(c, n)
   implicit none
   complex(kind=fd) c(*)
   integer n, i
   do i = 1, n
      c(i) = (0.d0, 0.d0)
   enddo
   return
end

! input
!	ik= wavenumber index (start at 1)
!	iklast= temporary wavenumber index to keep the bessel function in memory
!	nc= number of layer
!	nr= number of receiver
!	ns= number of source
!	nrs= number of radial distances
! output:
!	iklast updated
subroutine reflect0(ik, iklast, nc,nr,ns,nrs)
   use parameter
   use dimension1
   use dimension2

   implicit none

! Local
   integer          :: ir, ier, ik, iklast, ic, nrs,nc,nr,ns
   real(kind=fd)     :: fj0, arg
   real(kind=fd)     :: vy, gp, gs
   complex(kind=fd) :: ccv, cvp, cvs, cc

!     initialisations pour kr=0.

   if (ik .eq. 1) then
      call initc(u, 11*nr*ns)
   endif

! Compute bessel functions J0, J1, k1,k2,k3,k4,k5,k0
! iklast has been initialized at 0
   do ir = 1, nrs   ! loop over radial distances
      arg = rr(ir)*kr
      if (ik .gt. nkmax) then
         call ff01ad(fj0, vy, arg, 0)
         call ff02ad(fj1(ir), vy, arg, 0)
      else
         if (ik .gt. iklast) then
            call ff01ad(jj0(ik, ir), vy, arg, 0)
            call ff02ad(jj1(ik, ir), vy, arg, 0)
         endif
         fj0 = jj0(ik, ir)
         fj1(ir) = jj1(ik, ir)
      endif

      if (rr(ir) .ne. 0) then
         k0(ir) = kr*fj0
         k2(ir) = fj1(ir)/rr(ir)
         k1(ir) = k0(ir) - 2.*k2(ir)
         k4(ir) = k1(ir)/rr(ir)
         k3(ir) = -(2.*k4(ir) + kr2*fj1(ir))
      else
!                Lorsque rr=0. il faut utiliser
!                un developpement limite
         k1(ir) = 0.
         k2(ir) = kr/2.
         k3(ir) = 0.
         k4(ir) = 0.
      endif
      k5(ir) = k0(ir) - k2(ir)

   enddo

   if(ik.gt.iklast) then
!$OMP CRITICAL
     iklast=ik
!$OMP END CRITICAL
   endif
!               Calcul des nombres d onde verticaux

   do ic = 1, nc
      ccv = 1.+ai/(qp(ic) + qs(ic))
!     cvp=vp(ic)*ccv/(1.+.25/qp(ic)/qp(ic))/(1.-xlnf/qp(ic))
!     cvs=vs(ic)*ccv/(1.+.25/qs(ic)/qs(ic))/(1.-xlnf/qs(ic))
!__Futterman_____________________ !!!!!!!!!!!!!!!!!!!!!!!!
!
!      cvp=vp(ic)/(1-xlnf/pi/qp(ic)-ai/2./qp(ic))
!      cvs=vs(ic)/(1-xlnf/pi/qs(ic)-ai/2./qs(ic))
!
!__Konec Futtermana_______________!!!!!!!!!!!!!!!!!!!!!!!!

!__Const Q_Kjartansson___________ !!!!!!!!!!!!!!!!!!!!!!!!
!
      gp = atan(1./qp(ic))/pi
      cvp = vp(ic)*xlnf**gp/(1 - ai*tan(pi*gp/2))

      gs = atan(1./qs(ic))/pi
      cvs = vs(ic)*xlnf**gs/(1 - ai*tan(pi*gs/2))

      cka(ic) = omega/cvp
      ckb(ic) = omega/cvs
      ckb2(ic) = ckb(ic)*ckb(ic)
      cka2(ic) = cka(ic)*cka(ic)
      cc = cka2(ic) - kr2
      cnu(ic) = sqrt(cc)
      if (imag(cnu(ic)) .gt. 0.d0) cnu(ic) = -cnu(ic)
      cc = ckb2(ic) - kr2
      cgam(ic) = sqrt(cc)
      if (imag(cgam(ic)) .gt. 0.d0) cgam(ic) = -cgam(ic)
   end do
   do ic = 1, nc
      c2(ic) = kr*kr/ai/cnu(ic)
   enddo

   return
end
end module
