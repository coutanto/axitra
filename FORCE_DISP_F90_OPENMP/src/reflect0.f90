!******************************************************************************
!*                                                    AXITRA Version 3.0      *
!*                     SUBROUTINE REFLECT0                                    *
!*                                                                            *
!*    Calcul de coefficients dependant du nombre d'onde kr et de fonctions    *
!*    de Bessel.                                                              *
!*                                                                            *
!******************************************************************************

module reflect0m
contains
subroutine initc(c, n)
   implicit none
   complex(kind=8) c(*)
   integer n, i
   do i = 1, n
      c(i) = (0.d0, 0.d0)
   enddo
   return
end

subroutine reflect0(ik, nc, nr,ns,nrs)

   use dimension1
   use dimension2
   use parameter

   implicit none

   integer :: ir,ik,nc,nrs,nr,ns,ic
   real(kind=8) :: arg,gs,gp,vy
   complex(kind=8) :: cc,cvp,ccv,cvs

!     initialisations pour kr=0.

   if (ik .eq. 1) then
      call initc(u, 5*nr*ns)
   endif

!     Calcul des fonctions de Bessel J0 et J1, k1,k2,k3,k4,k5,k0

   do ir = 1, nrs

      arg = rr(ir)*kr

      call ff01ad(j0(ir), vy, arg, 0)
      call ff02ad(j1(ir), vy, arg, 0)

      if (rr(ir) .ne. 0) then
         k0(ir) = kr*j0(ir)
         k2(ir) = j1(ir)/rr(ir)
         k1(ir) = k0(ir) - 2.*k2(ir)
         k4(ir) = k1(ir)/rr(ir)
         k3(ir) = -(2.*k4(ir) + kr2*j1(ir))
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

!               Calcul des nombres d'onde verticaux

   do ic = 1, nc
#define OLDEF
#ifdef OLDEF !old attenuation definition
      ccv = 1.+ai/(qp(ic) + qs(ic))
      cvp = vp(ic)*ccv/(1.+.25/qp(ic)/qp(ic))/(1.-xlnf/qp(ic))
      cvs = vs(ic)*ccv/(1.+.25/qs(ic)/qs(ic))/(1.-xlnf/qs(ic))
#else ! new definition implemented by Manu Chaljub
      gp = atan(1./qp(ic))/pi
      cvp = vp(ic)*xlnf**gp/(1 - ai*tan(pi*gp/2))

      gs = atan(1./qs(ic))/pi
      cvs = vs(ic)*xlnf**gs/(1 - ai*tan(pi*gs/2))
#endif
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
   enddo
   do ic = 1, nc
      c2(ic) = kr*kr/ai/cnu(ic)
!write(0,*) kr,cnu(1),cvp
   enddo

   return
end
end module
