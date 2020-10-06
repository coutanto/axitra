!******************************************************************************
!*                     AXITRA Force Version
!                         SUBROUTINE REFLECT1
!
!               Calcul des coefficients de reflexion/transmission
!                   Matrice de Reflexion/Transmission et Dephasage
!       (Les coefficients de reflexion/transmission utilisent les memes
!            termes intermediaires que Aki-Richards p149, MAIS :
!            Aki utilise la convention inverse pour la TF (exp(-iwt)),
!        et travaille avec le parametre de rai et les angles d'incidences)
!
!      Le potentiel PSI utilise pour l'onde SV est defini par :
!                  u = rot ( rot (PSI) )
!      i.e. un terme de derivation supplementaire par rapport a la convention
!      habituelle : u= rot (PSI)
!
!   On deduit les coefficients de REF/TRANS de ceux definis par la convention
!      classique en divisant le potentiel PSI par 1./ai/kr = coef
!
!       Ordre de stockage :
!                              pp=(1,1)   sp=(1,2)
!                              ps=(2,1)   ss=(2,2)                            
!******************************************************************************
module reflect1m

contains
subroutine reflect1(freesurface, nc, uflow)

   use parameter
   use dimension1
   use dimension2

   implicit none

   logical :: freesurface,uflow
   integer :: nc

   integer :: ic,ic1
   complex(kind=8) :: coef, cf1, cf2, cf3, cdd, ca, cb, cb1, cb2, ca1d, arg
   complex(kind=8) :: ca2d, cc, cd, ce, cf, cg, ch, cdeph, cs1, cs2, cdelt
   real(kind=8) :: aki
! Coefficient pour la convention sur PSI (coef) et sur la TF (aki=-1.)
   coef = 1./ai
   aki = -1.

!               CONDITIONS AUX LIMITES a la profondeur z=0, coefficients
!               de reflexion/transmission
!      2 possibilites : 1) surface libre
!                       2) 1/2 espace sup. infini
   ru=0.d0;rd=0.d0;tu=0.d0;td=0.d0
!A1                    SURFACE LIBRE
   if (freesurface) then

      cf1 = (ckb2(1) - 2.*kr2)
      cf2 = cf1*cf1
      cf3 = 4.*cnu(1)*kr2*cgam(1)
      cdd = cf2 + cf3

      ru(1, 1, 1) = (-cf2 + cf3)/cdd
      ru(1, 2, 1) = 4.*cnu(1)*cf1/cdd*coef*aki
      ru(1, 2, 2) = (cf2 - cf3)/cdd*aki
      ru(1, 1, 2) = 4.*kr2*cgam(1)*cf1/cdd/coef
      tu(1, 1, 1) = 0.
      tu(1, 1, 2) = 0.
      tu(1, 2, 1) = 0.
      tu(1, 2, 2) = 0.
      rush(1) = 1.
      tush(1) = 0.
!A2
   else
!B1                   1/2 ESPACE SUP. INFINI
      ru(1, 1, 1) = 0.
      ru(1, 2, 1) = 0.
      ru(1, 2, 2) = 0.
      ru(1, 1, 2) = 0.
      tu(1, 1, 1) = 1.
      tu(1, 1, 2) = 0.
      tu(1, 2, 1) = 0.
      tu(1, 2, 2) = 1.
      rush(1) = 0.
      tush(1) = 1.
!B2
   endif
!               Coefficients aux interfaces entre couches

   do ic = 2, nc

      ic1 = ic -1
      cb1 = kr2/ckb2(ic1)
      cb2 = kr2/ckb2(ic)
      ca1d = rho(ic1)*(1.-2.*cb1)
      ca2d = rho(ic)*(1.-2.*cb2)
      ca = ca2d-ca1d
      cb = ca2d+2.*rho(ic1)*cb1
      cc = ca1d+2.*rho(ic)*cb2
      cd = 2.*(rho(ic)/ckb2(ic) - rho(ic1)/ckb2(ic1))
      ce = cb*cnu(ic1) + cc*cnu(ic)
      cf = cb*cgam(ic1) + cc*cgam(ic)
      cg = ca - cd*cnu(ic1)*cgam(ic)
      ch = ca - cd*cnu(ic)*cgam(ic1)
      cdd = ce*cf + cg*ch*kr2

      rd(ic, 1, 1) = (cf*(cb*cnu(ic1) - cc*cnu(ic)) - ch*kr2*(ca + cd*cnu(ic1)*cgam(ic)))/cdd
      rd(ic, 1, 2) = -2.*kr2*cgam(ic1)*(ca*cb + cc*cd*cnu(ic)*cgam(ic))/cdd/coef*aki
      rd(ic, 2, 2) = -(ce*(cb*cgam(ic1) - cc*cgam(ic)) - cg*kr2*(ca + cd*cnu(ic)*cgam(ic1)))/cdd*aki
      rd(ic, 2, 1) = -2.*cnu(ic1)*(ca*cb + cc*cd*cnu(ic)*cgam(ic))/cdd*coef
      td(ic, 1, 1) = 2.*rho(ic1)*cnu(ic1)*cf/cdd
      td(ic, 1, 2) = -2.*rho(ic1)*cgam(ic1)*cg*kr2/cdd/coef*aki
      td(ic, 2, 2) = 2.*rho(ic1)*cgam(ic1)*ce/cdd
      td(ic, 2, 1) = 2.*rho(ic1)*cnu(ic1)*ch/cdd*coef*aki

      ru(ic, 1, 1) = -(cf*(cb*cnu(ic1) - cc*cnu(ic)) + cg*kr2*(ca + cd*cnu(ic)*cgam(ic1)))/cdd
      ru(ic, 1, 2) = 2.*kr2*cgam(ic)*(ca*cc + cb*cd*cnu(ic1)*cgam(ic1))/cdd/coef
      ru(ic, 2, 2) = (ce*(cb*cgam(ic1) - cc*cgam(ic)) + ch*kr2*(ca + cd*cnu(ic1)*cgam(ic)))/cdd*aki
      ru(ic, 2, 1) = 2.*cnu(ic)*(ca*cc + cb*cd*cnu(ic1)*cgam(ic1))/cdd*coef*aki
      tu(ic, 1, 1) = 2.*rho(ic)*cnu(ic)*cf/cdd
      tu(ic, 1, 2) = 2.*rho(ic)*cgam(ic)*ch*kr2/cdd/coef
      tu(ic, 2, 2) = 2.*rho(ic)*cgam(ic)*ce/cdd
      tu(ic, 2, 1) = -2.*rho(ic)*cnu(ic)*cg/cdd*coef

!   Modification pour calculateur a faible dynamique [1.e-300; 1.e+300]
      arg = -ai*cnu(ic1)*hc(ic1)
      if (dreal(arg) .lt. explim) then !exp(-inf)*exp(i*theta)=0.
         uflow=.true.
         me1(ic1)=0.d0
      else
         me1(ic1) = exp(arg)
      endif
      arg   = -ai*cgam(ic1)*hc(ic1)
      if (dreal(arg) .lt. explim) then
        uflow=.true.
        me2(ic1) = 0.d0
      else
        me2(ic1) = exp(arg)
      endif

      cs1 = rho(ic1)/ckb2(ic1)*cgam(ic1)
      cs2 = rho(ic)/ckb2(ic)*cgam(ic)
      cdelt = cs1 + cs2

      rush(ic) = (cs2 - cs1)/cdelt
      rdsh(ic) = -rush(ic)
      tush(ic) = 2.*cs2/cdelt
      tdsh(ic) = 2.*cs1/cdelt

   enddo

   return
end
end module
