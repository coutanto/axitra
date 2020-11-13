!***********************************************************
!                     AXITRA Moment Version
!*
!*              SUBROUTINE REFLECT5
!*
!*        Calcul des deplacements avec diverses rotations
!*        et recombinaisons. Passage aux sources du tenseur
!*        des moments sismiques (M1 a M6)
!*        Multiplication par les termes frequentiel
!*        et angulaire :
!*                     u=u*C3(theta)*CFF(omega )           
!***********************************************************
module reflect5m
contains
subroutine reflect5(jf,nr,ns)

   use dimension1
   use dimension2
   use parameter
   implicit none

   integer base, is, ir, it, jf, nr,ns, i
   real(kind=fd)     :: cor, cor2, co2r, sir, sir2, si2r
   complex(kind=fd) :: urxx, utxx, uzxx, urxy, utxy, uzxy, urxz, utxz, uzxz
   complex(kind=fd) :: uryy, utyy, uzyy, uryz, utyz, uzyz, urzz, utzz, uzzz
   complex(kind=fd) :: ux(6), uy(6), uz(6)
   complex(kind=fd),allocatable :: buffer(:,:,:)

   allocate(buffer(18,nr,ns))

   do is = 1, ns
   do ir = 1, nr
   do it = 1, 11
      u(ir, is, it) = u(ir, is, it)*a1*cff(is)
   enddo
!+++++++++++++
!        Deplacement dus aux sources Mxx,Mxy,Mxz,Myy,Myz,Mzz avec
!        convention de signe inverse pour le deplacement vertical
!        (positif vers le haut)
!+++++++++++++

   cor = cosr(ir, is)
   sir = sinr(ir, is)
   cor2 = cor*cor
   sir2 = sir*sir
   co2r = cor2 - sir2
   si2r = 2.*cor*sir

!        Mxx

   urxx = -cor2*u(ir, is, 1) - u(ir, is, 2) - co2r*u(ir, is, 5)
   utxx = si2r*(u(ir, is, 2) + u(ir, is, 6)/2.)
   uzxx = cor2*u(ir, is, 3) + u(ir, is, 4)

!        Mxy+Myx

   urxy = -si2r*(u(ir, is, 1) + 2.*u(ir, is, 5))
   utxy = -co2r*(2.*u(ir, is, 2) + u(ir, is, 6))
   uzxy = si2r*u(ir, is, 3)

!        Mxz+Mzx

   urxz = -cor*u(ir, is, 7)
   utxz = sir*u(ir, is, 8)
   uzxz = cor*u(ir, is, 9)

!        Myy

   uryy = -sir2*u(ir, is, 1) - u(ir, is, 2) + co2r*u(ir, is, 5)
   utyy = -si2r*(u(ir, is, 2) + u(ir, is, 6)/2.)
   uzyy = sir2*u(ir, is, 3) + u(ir, is, 4)

!        Myz+Mzy

   uryz = -sir*u(ir, is, 7)
   utyz = -cor*u(ir, is, 8)
   uzyz = sir*u(ir, is, 9)

!        Mzz

   urzz = -u(ir, is, 10)
   utzz = 0.
   uzzz = -u(ir, is, 11)

!+++++++++++
!     Passage aux sources bis, 5 dislocations elementaires et
!     une source isotrope + rotation des composantes pour passer de
!     radial/tangentiel a Ox/Oy
!+++++++++++
!        M1 = (Mxy + Myx)
   ux(1) = urxy*cor - utxy*sir
   uy(1) = urxy*sir + utxy*cor
   uz(1) = uzxy

!        M2 = (Mxz + Mzx)
   ux(2) = urxz*cor - utxz*sir
   uy(2) = urxz*sir + utxz*cor
   uz(2) = uzxz

!        M3 = -(Myz + Mzy)
   ux(3) = -uryz*cor + utyz*sir
   uy(3) = -uryz*sir - utyz*cor
   uz(3) = -uzyz

!        M4 = -Mxx + Mzz
   ux(4) = (-urxx + urzz)*cor - (-utxx + utzz)*sir
   uy(4) = (-urxx + urzz)*sir + (-utxx + utzz)*cor
   uz(4) = -uzxx + uzzz

!        M5 = -Myy + Mzz
   ux(5) = (-uryy + urzz)*cor - (-utyy + utzz)*sir
   uy(5) = (-uryy + urzz)*sir + (-utyy + utzz)*cor
   uz(5) = -uzyy + uzzz

!        M6 = Mxx + Myy + Mzz
   ux(6) = (urxx + uryy + urzz)*cor - (utxx + utyy + utzz)*sir
   uy(6) = (urxx + uryy + urzz)*sir + (utxx + utyy + utzz)*cor
   uz(6) = uzxx + uzyy + uzzz

!   base = (ir - 1)*3 + (is - 1)*3*nr + (jf - 1)*3*nr*ns
!   write (out2, rec=base + 1) (ux(it), it=1, 6)
!   write (out2, rec=base + 2) (uy(it), it=1, 6)
!   write (out2, rec=base + 3) (uz(it), it=1, 6)

   do i=1,6
     buffer(i,ir,is)=ux(i)
     buffer(i+6,ir,is)=uy(i)
     buffer(i+12,ir,is)=uz(i)
   enddo

   enddo ! receiver loop
   enddo ! source loop

   write(out2,rec=jf) buffer
   deallocate(buffer)
   return
end
end module
