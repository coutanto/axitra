!*******************************************************************************
!*                     AXITRA Force Version
!*
!*                         SUBROUTINE REFLECT2
!*
!*       Calcul des matrices de reflectivite mt,mb ou nt,nb dans chaque couche
!*         (rapport des potentiels montant/descendant ou descendant/montant)
!*       Le suffixe t ou b precise si la matrice est donnee au sommet (top)
!*       ou au bas (bottom) de la couche.
!*       fup et fdo sont des matrices intermediaires utilisees dans le calcul
!*       des potentiels.
!*       Ordre de stockage :
!*                    pp=(1,1)   sp=(1,2)
!*                    ps=(2,1)   ss=(2,2)
!*******************************************************************************
module reflect2m

contains
subroutine reflect2(nc)

   use dimension1
   use dimension2
   use parameter
   implicit none

   complex(kind=8) :: nb(2, 2), mb(2, 2), nbsh, mbsh
   complex(kind=8) :: ca1, ca2, ca3, ca4, cb1, cb2, cb3, cb4, cc1, cc2, cc3, cc4
   complex(kind=8) :: cadet, cash, cbsh, ccsh
   integer          :: ic, ic1, nc
!
!                        Calcul pour les couches au dessus de la source
!

   nt(1, 1, 1) = ru(1, 1, 1)
   nt(1, 1, 2) = ru(1, 1, 2)
   nt(1, 2, 1) = ru(1, 2, 1)
   nt(1, 2, 2) = ru(1, 2, 2)
   ntsh(1) = rush(1)

   do ic = 1, nc - 1

      nb(1, 1) = me1(ic)*me1(ic)*nt(ic, 1, 1)
      nb(1, 2) = me1(ic)*me2(ic)*nt(ic, 1, 2)
      nb(2, 1) = me2(ic)*me1(ic)*nt(ic, 2, 1)
      nb(2, 2) = me2(ic)*me2(ic)*nt(ic, 2, 2)
      nbsh = me2(ic)*me2(ic)*ntsh(ic)

      ca1 = 1.-(rd(ic + 1, 1, 1)*nb(1, 1) + rd(ic + 1, 1, 2)*nb(2, 1))
      ca2 = -(rd(ic + 1, 1, 1)*nb(1, 2) + rd(ic + 1, 1, 2)*nb(2, 2))
      ca3 = -(rd(ic + 1, 2, 1)*nb(1, 1) + rd(ic + 1, 2, 2)*nb(2, 1))
      ca4 = 1.-(rd(ic + 1, 2, 1)*nb(1, 2) + rd(ic + 1, 2, 2)*nb(2, 2))
      cadet = ca1*ca4 - ca2*ca3
      cash = 1./(1.-rdsh(ic + 1)*nbsh)

      cb1 = td(ic + 1, 1, 1)*nb(1, 1) + td(ic + 1, 1, 2)*nb(2, 1)
      cb2 = td(ic + 1, 1, 1)*nb(1, 2) + td(ic + 1, 1, 2)*nb(2, 2)
      cb3 = td(ic + 1, 2, 1)*nb(1, 1) + td(ic + 1, 2, 2)*nb(2, 1)
      cb4 = td(ic + 1, 2, 1)*nb(1, 2) + td(ic + 1, 2, 2)*nb(2, 2)
      cbsh = tdsh(ic + 1)*nbsh

      cc1 = (ca4*tu(ic + 1, 1, 1) - ca2*tu(ic + 1, 2, 1))/cadet
      cc2 = (ca4*tu(ic + 1, 1, 2) - ca2*tu(ic + 1, 2, 2))/cadet
      cc3 = (-ca3*tu(ic + 1, 1, 1) + ca1*tu(ic + 1, 2, 1))/cadet
      cc4 = (-ca3*tu(ic + 1, 1, 2) + ca1*tu(ic + 1, 2, 2))/cadet
      ccsh = cash*tush(ic + 1)

      nt(ic + 1, 1, 1) = ru(ic + 1, 1, 1) + cb1*cc1 + cb2*cc3
      nt(ic + 1, 1, 2) = ru(ic + 1, 1, 2) + cb1*cc2 + cb2*cc4
      nt(ic + 1, 2, 1) = ru(ic + 1, 2, 1) + cb3*cc1 + cb4*cc3
      nt(ic + 1, 2, 2) = ru(ic + 1, 2, 2) + cb3*cc2 + cb4*cc4
      ntsh(ic + 1) = rush(ic + 1) + cbsh*ccsh

      fup(ic, 1, 1) = cc1*me1(ic)
      fup(ic, 1, 2) = cc2*me1(ic)
      fup(ic, 2, 1) = cc3*me2(ic)
      fup(ic, 2, 2) = cc4*me2(ic)
      fupsh(ic) = ccsh*me2(ic)

   enddo
!
!                        Calcul pour les couches au dessous de la source
!

   mt(nc, 1, 1) = 0.
   mt(nc, 1, 2) = 0.
   mt(nc, 2, 1) = 0.
   mt(nc, 2, 2) = 0.
   mtsh(nc) = 0.

   do ic = nc - 1, 1, -1

      ca1 = 1.-(ru(ic + 1, 1, 1)*mt(ic + 1, 1, 1) + ru(ic + 1, 1, 2)*mt(ic + 1, 2, 1))
      ca2 = -(ru(ic + 1, 1, 1)*mt(ic + 1, 1, 2) + ru(ic + 1, 1, 2)*mt(ic + 1, 2, 2))
      ca3 = -(ru(ic + 1, 2, 1)*mt(ic + 1, 1, 1) + ru(ic + 1, 2, 2)*mt(ic + 1, 2, 1))
      ca4 = 1.-(ru(ic + 1, 2, 1)*mt(ic + 1, 1, 2) + ru(ic + 1, 2, 2)*mt(ic + 1, 2, 2))
      cadet = ca1*ca4 - ca2*ca3
      cash = 1./(1.-rush(ic + 1)*mtsh(ic + 1))

      cb1 = tu(ic + 1, 1, 1)*mt(ic + 1, 1, 1) + tu(ic + 1, 1, 2)*mt(ic + 1, 2, 1)
      cb2 = tu(ic + 1, 1, 1)*mt(ic + 1, 1, 2) + tu(ic + 1, 1, 2)*mt(ic + 1, 2, 2)
      cb3 = tu(ic + 1, 2, 1)*mt(ic + 1, 1, 1) + tu(ic + 1, 2, 2)*mt(ic + 1, 2, 1)
      cb4 = tu(ic + 1, 2, 1)*mt(ic + 1, 1, 2) + tu(ic + 1, 2, 2)*mt(ic + 1, 2, 2)
      cbsh = tush(ic + 1)*mtsh(ic + 1)

      cc1 = (ca4*td(ic + 1, 1, 1) - ca2*td(ic + 1, 2, 1))/cadet
      cc2 = (ca4*td(ic + 1, 1, 2) - ca2*td(ic + 1, 2, 2))/cadet
      cc3 = (-ca3*td(ic + 1, 1, 1) + ca1*td(ic + 1, 2, 1))/cadet
      cc4 = (-ca3*td(ic + 1, 1, 2) + ca1*td(ic + 1, 2, 2))/cadet
      ccsh = cash*tdsh(ic + 1)

      mb(1, 1) = rd(ic + 1, 1, 1) + cb1*cc1 + cb2*cc3
      mb(1, 2) = rd(ic + 1, 1, 2) + cb1*cc2 + cb2*cc4
      mb(2, 1) = rd(ic + 1, 2, 1) + cb3*cc1 + cb4*cc3
      mb(2, 2) = rd(ic + 1, 2, 2) + cb3*cc2 + cb4*cc4
      mbsh = rdsh(ic + 1) + cbsh*ccsh

      mt(ic, 1, 1) = me1(ic)*me1(ic)*mb(1, 1)
      mt(ic, 1, 2) = me1(ic)*me2(ic)*mb(1, 2)
      mt(ic, 2, 1) = me2(ic)*me1(ic)*mb(2, 1)
      mt(ic, 2, 2) = me2(ic)*me2(ic)*mb(2, 2)
      mtsh(ic) = me2(ic)*me2(ic)*mbsh

      fdo(ic + 1, 1, 1) = cc1*me1(ic)
      fdo(ic + 1, 1, 2) = cc2*me2(ic)
      fdo(ic + 1, 2, 1) = cc3*me1(ic)
      fdo(ic + 1, 2, 2) = cc4*me2(ic)
      fdosh(ic + 1) = ccsh*me2(ic)

   enddo

   return
end
end module
