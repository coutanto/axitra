!***********************************************************
!*                                      AXITRA Version 3.0 *
!*              SUBROUTINE REFLECT5V3                      *
!*                                                         *
!*        Calcul des deplacements avec diverses rotations  *
!*        et recombinaisons.                               *
!*        Multiplication par les termes frequentiel        *
!*        et angulaire :                                   *
!*                     u=u*C3(theta)*CFF(omega )           *
!***********************************************************
module reflect5m

contains
subroutine reflect5(jf,nr,ns)

   use dimension1
   use dimension2
   use parameter

   implicit none

   integer          :: it,jf,base,is,nr,ns,ir,i
   real(kind=fd)     :: cor, sir, cof
   complex(kind=fd) :: ur(3), ut(3), uz(3),a, ux(3),uy(3)
   complex(kind=fd),allocatable :: buffer(:,:,:)

   allocate(buffer(9,nr,ns))

   do is = 1, ns
      do ir = 1, nr
	 a=a1*cff(is)
         u(ir, is, 1) = u(ir, is, 1)*a
         u(ir, is, 2) = u(ir, is, 2)*a
         u(ir, is, 3) = u(ir, is, 3)*a
         u(ir, is, 4) = u(ir, is, 4)*a
         u(ir, is, 5) = u(ir, is, 5)*a

!+++++++++++++
!        Displacements due to forces Fx, Fy, Fz
!        in radial, transverse, vertical
!+++++++++++++
#define MV08052018   ! correction proposed by Martin Vallée
#ifdef MV08052018
         cof=-1_8
#else
         cof=1._8
#endif

         cor = cosr(ir, is)
         sir = sinr(ir, is)

!        Fx
         ur(1) = cor*u(ir, is, 1)
         ut(1) = -sir*u(ir, is, 2)
         uz(1) = cof*cor*u(ir, is, 3)

!        Fy
         ur(2) = sir*u(ir, is, 1)
         ut(2) = cor*u(ir, is, 2)
         uz(2) = cof*sir*u(ir, is, 3)

!        Fz
         ur(3) = cof*u(ir, is, 4)
         ut(3) = 0.
         uz(3) = cof*cof*u(ir, is, 5)
!
!  and in X,Y, Z
!
         do i=1,3
           ux(i) = ur(i)*cor - ut(i)*sir
           uy(i) = ur(i)*sir + ut(i)*cor
         enddo

!		 base = (ir - 1)*3 + (is - 1)*3*nr + (jf - 1)*3*nr*ns
!         write (out2, rec=base + 1) (ux(it), it=1, 3)
!         write (out2, rec=base + 2) (uy(it), it=1, 3)
!         write (out2, rec=base + 3) (uz(it), it=1, 3)
         buffer(1,ir,is)=ux(1)
         buffer(2,ir,is)=ux(2)
         buffer(3,ir,is)=ux(3)
         buffer(4,ir,is)=uy(1)
         buffer(5,ir,is)=uy(2)
         buffer(6,ir,is)=uy(3)
         buffer(7,ir,is)=uz(1)
         buffer(8,ir,is)=uz(2)
         buffer(9,ir,is)=uz(3)
      enddo
   enddo
   write(out2,rec=jf) buffer
   deallocate(buffer)
   return
end
end module
