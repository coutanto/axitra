!**************************************************************
!                     AXITRA Force Version
!*
!                  SUBROUTINE REFLECT3
!
!  Calcul des potentiels dus a 6 sources elementaires, au
!  sommet de la couche source ISC. Les 6 sources elementaires
!  sont 6 sources de potentiels, 2 PHI, 2 PSI, 2 KHI. Ces
!  sources different par l'existence d'un terme sign(z-z0).
!**************************************************************
module reflect3m

contains
subroutine reflect3(ncs)

   use parameter
   use dimension1
   use dimension2

   implicit none

   integer is0, is1, ic, is2, ncs
   real(kind=8)     :: zsc
   complex(kind=8) :: rud1, rud2, rud3, rud4, rdu1, rdu2, rdu3, rdu4, tud1
   complex(kind=8) :: tud2, tud3, tud4, tdu1, tdu2, tdu3, tdu4, tsh, egam
   complex(kind=8) :: enu, rupsh, rdosh, arg, cdet, cu1sh, cd1sh, cu2sh, cd2sh

   complex(kind=8) :: cu1(2), cd1(2), cu2(2), cd2(2), cu3(2), cd3(2)
   complex(kind=8) :: cu4(2), cd4(2), rup(2, 2), rdo(2, 2)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Matrice de passage des vecteurs potentiels sources su0, sd0
!  aux vecteurs potentiel de la couche source isc: su, sd
!
!                    [tud] et [tdu]
!
!                 ------------------------
!
!  su(,) : potentiel montant au sommet de la couche
!  sd(,) : potentiel descendant au sommet de la couche
!
!     (*,) : type de source (1 a 5)
!     (,*) : type de potentiel PHI ou PSI=KHI (1, 2)
!
!                 ------------------------
!  Les vecteurs potentiels su() et sd() sont obtenus a
!  partir des potentiels des sources su0(), sd0() au
!  sommet de la couche source par :
!
!     su = 1/[1 - rup*rdo] . (su0 + [rup].sd0)
!
!     sd = 1/[1 - rdo*rup] . (sd0 + [rdo].su0)
!
!
!     ou les matrices rup et rdo sont donnees par les
! matrices reflectivite du sommet de la couche source isc :
!
!                [rup] = [mt(isc)]
!                [rdo] = [nt(isc)]
!
!        [rdo] = matrice reflectivite DOWN
!           (potentiel descendant/potentiel montant) due
!     a l'empilement de couches situe au dessus de la source
!
!        [rup] = matrice reflectivite UP
!           (potentiel montant/potentiel descendant) due
!     a l'empilement de couches situe au dessous de la source
!
!        [rud] = [rup] * [rdo]
!
!        [rdu] = [rdo] * [rup]
!
!     On pose [tud] = 1/[1 - rup*rdo]
!             [tdu] = 1/[1 - rdo*rup]
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do is1 = 1, ncs !boucle sur les couches sources
      ic = isc(is1)
      rdo(1, 1) = nt(ic, 1, 1)
      rdo(1, 2) = nt(ic, 1, 2)
      rdo(2, 1) = nt(ic, 2, 1)
      rdo(2, 2) = nt(ic, 2, 2)
      rdosh = ntsh(ic)
      rup(1, 1) = mt(ic, 1, 1)
      rup(1, 2) = mt(ic, 1, 2)
      rup(2, 1) = mt(ic, 2, 1)
      rup(2, 2) = mt(ic, 2, 2)
      rupsh = mtsh(ic)

      rud1 = rup(1, 1)*rdo(1, 1) + rup(1, 2)*rdo(2, 1)
      rud2 = rup(1, 1)*rdo(1, 2) + rup(1, 2)*rdo(2, 2)
      rud3 = rup(2, 1)*rdo(1, 1) + rup(2, 2)*rdo(2, 1)
      rud4 = rup(2, 1)*rdo(1, 2) + rup(2, 2)*rdo(2, 2)

      rdu1 = rdo(1, 1)*rup(1, 1) + rdo(1, 2)*rup(2, 1)
      rdu2 = rdo(1, 1)*rup(1, 2) + rdo(1, 2)*rup(2, 2)
      rdu3 = rdo(2, 1)*rup(1, 1) + rdo(2, 2)*rup(2, 1)
      rdu4 = rdo(2, 1)*rup(1, 2) + rdo(2, 2)*rup(2, 2)

      cdet = (1.-rud1)*(1.-rud4) - rud2*rud3

      tud1 = (1.-rud4)/cdet
      tud2 = rud2/cdet
      tud3 = rud3/cdet
      tud4 = (1.-rud1)/cdet

      cdet = (1.-rdu1)*(1.-rdu4) - rdu2*rdu3

      tdu1 = (1.-rdu4)/cdet
      tdu2 = rdu2/cdet
      tdu3 = rdu3/cdet
      tdu4 = (1.-rdu1)/cdet

      tsh = 1./(1.-rupsh*rdosh)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Vecteurs potentiel source pour 4 sources elementaires :
!        (dephasage calcule / sommet de la couche)
!
!              cui = su0 + [rup].sd0  (i=1,4)
!              cdi = sd0 + [rdo].su0
!
!  et potentiel KHI couche source pour 2 sources elementaires :
!
!                      cuish = su0sh + rupsh*sd0sh (i=1,2)
!                      cdish = sd0sh + rdosh*su0sh
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do is2 = 1, nzs(is1) !boucle sur les prof sources
         zsc = zs(izss(1, is2, is1))
         is0 = izss(1, is2, is1)
         arg = -ai*cgam(ic)*zsc
         if (real(arg) .lt. explim) arg = cmplx(explim, imag(arg))
         egam = exp(arg)
         arg = -ai*cnu(ic)*zsc
         if (real(arg) .lt. explim) arg = cmplx(explim, imag(arg))
         enu = exp(arg)

!                        Source PHI
         cu1(1) = enu + rup(1, 1)/enu
         cu1(2) = rup(2, 1)/enu
         cd1(1) = 1./enu + rdo(1, 1)*enu
         cd1(2) = rdo(2, 1)*enu
!                        Source PHI*sign(z-z0)
         cu2(1) = -enu + rup(1, 1)/enu
         cu2(2) = rup(2, 1)/enu
         cd2(1) = 1./enu - rdo(1, 1)*enu
         cd2(2) = -rdo(2, 1)*enu
!                        Source PSI
         cu3(1) = rup(1, 2)/egam
         cu3(2) = egam + rup(2, 2)/egam
         cd3(1) = rdo(1, 2)*egam
         cd3(2) = 1./egam + rdo(2, 2)*egam
!                        Source PSI*sign(z-z0)
         cu4(1) = rup(1, 2)/egam
         cu4(2) = -egam + rup(2, 2)/egam
         cd4(1) = -rdo(1, 2)*egam
         cd4(2) = 1./egam - rdo(2, 2)*egam
!                        Source KHI
         cu1sh = egam + rupsh/egam
         cd1sh = 1./egam + rdosh*egam
!                        Source KHI*sign(z-z0)
         cu2sh = -egam + rupsh/egam
         cd2sh = 1./egam - rdosh*egam

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Potentiels PHI, PSI et KHI, montant et descendant, dans la couche
!   source (dephasage / sommet) pour les 6 sources elementaires.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                        Source PHI
         su1(is0, 1) = tud1*cu1(1) + tud2*cu1(2)
         su1(is0, 2) = tud3*cu1(1) + tud4*cu1(2)
         sd1(is0, 1) = tdu1*cd1(1) + tdu2*cd1(2)
         sd1(is0, 2) = tdu3*cd1(1) + tdu4*cd1(2)
!                        Source PHI*sign(z-z0)
         su2(is0, 1) = tud1*cu2(1) + tud2*cu2(2)
         su2(is0, 2) = tud3*cu2(1) + tud4*cu2(2)
         sd2(is0, 1) = tdu1*cd2(1) + tdu2*cd2(2)
         sd2(is0, 2) = tdu3*cd2(1) + tdu4*cd2(2)
!                        Source PSI
         su3(is0, 1) = tud1*cu3(1) + tud2*cu3(2)
         su3(is0, 2) = tud3*cu3(1) + tud4*cu3(2)
         sd3(is0, 1) = tdu1*cd3(1) + tdu2*cd3(2)
         sd3(is0, 2) = tdu3*cd3(1) + tdu4*cd3(2)
!                        Source PSI*sign(z-z0)
         su4(is0, 1) = tud1*cu4(1) + tud2*cu4(2)
         su4(is0, 2) = tud3*cu4(1) + tud4*cu4(2)
         sd4(is0, 1) = tdu1*cd4(1) + tdu2*cd4(2)
         sd4(is0, 2) = tdu3*cd4(1) + tdu4*cd4(2)
!                        Source KHI
         su1sh(is0) = tsh*cu1sh
         sd1sh(is0) = tsh*cd1sh
!                        Source KHI*sign(z-z0)
         su2sh(is0) = tsh*cu2sh
         sd2sh(is0) = tsh*cd2sh
      enddo
   enddo
   return
end
end module
