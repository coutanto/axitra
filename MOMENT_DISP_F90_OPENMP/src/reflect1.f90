!******************************************************************************
!                     AXITRA Moment Version
!                                                                             *
!                         SUBROUTINE REFLECT1                                 *
!                                                                             *
!               Calcul des coefficients de reflexion/transmission             *
!                   Matrice de Reflexion/Transmission et Dephasage             *
!       (Les coefficients de reflexion/transmission utilisent les memes       *
!            termes intermediaires que Aki-Richards p149, MAIS :              *
!            Aki utilise la convention inverse pour la TF (exp(-iwt)),        *
!        et travaille avec le parametre de rai et les angles d incidences)    *
!                                                                             *
!      Le potentiel PSI utilise pour l onde SV est defini par :               *
!                  u = rot ( rot (PSI) )                                      *
!      i.e. un terme de derivation supplementaire par rapport a la convention *
!      habituelle : u= rot (PSI)                                              *
!                                                                             *
!   On deduit les coefficients de REF/TRANS de ceux definis par la convention *
!      classique en divisant le potentiel PSI par 1./ai/kr = coef             *
!                                                                             *
!       Ordre de stockage :                                                   *
!                              pp=(1,1)   sp=(1,2)                            *
!                              ps=(2,1)   ss=(2,2)                            *
!******************************************************************************
module reflect1m
contains
subroutine reflect1(freeSurface,nc,uflow)

   use parameter
   use dimension1
   use dimension2

   implicit none
   logical :: freeSurface

   integer         :: ic, ic1,nc
   real(kind=8)    :: aki
   complex(kind=8) :: coef, cf1, cf2, cf3, cdd, cnurho, cgarho, cnugam, cb1, cb2, ckb2i1, ckb2i2
   complex(kind=8) :: ca1d, ca2d, ca, cb, cc, cd, ce, cf, cg, ch, cgkr2, chkr2, cdd2, cs1, cs2, cdelt
   complex(kind=8) :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7, ctmp8, ctmp9, ctmp10, arg
   logical         :: uflow

! Coefficient pour la convention sur PSI (coef) et sur la TF (aki=-1.)
   coef = 1./ai
   aki = -1.

!               CONDITIONS AUX LIMITES a la profondeur z=0, coefficients
!               de reflexion/transmission
!      2 possibilites : 1) surface libre
!                       2) 1/2 espace sup. infini 
   ru=0.d0;rd=0.d0;tu=0.d0;td=0.d0
!A1                    SURFACE LIBRE
   if (freeSurface) then
      cf1 = (ckb2(1) - 2.*kr2)
      cf2 = cf1*cf1
      cf3 = 4.d0*cnu(1)*kr2*cgam(1)
      cdd = 1.d0/(cf2 + cf3)

      ru(1, 1, 1) = (-cf2 + cf3)*cdd
      ru(1, 2, 1) = 4.d0*cnu(1)*cf1*cdd*coef*aki
      ru(1, 2, 2) = (cf2 - cf3)*cdd*aki
      ru(1, 1, 2) = 4.d0*kr2*cgam(1)*cf1*cdd*ai
      tu(1, 1, 1) = 0.d0
      tu(1, 1, 2) = 0.d0
      tu(1, 2, 1) = 0.d0
      tu(1, 2, 2) = 0.d0
      rush(1) = 1.d0
      tush(1) = 0.d0
   else

!B1                   1/2 ESPACE SUP. INFINI
      ru(1, 1, 1) = 0.d0
      ru(1, 2, 1) = 0.d0
      ru(1, 2, 2) = 0.d0
      ru(1, 1, 2) = 0.d0
      tu(1, 1, 1) = 1.d0
      tu(1, 1, 2) = 0.d0
      tu(1, 2, 1) = 0.d0
      tu(1, 2, 2) = 1.d0
      rush(1) = 0.d0
      tush(1) = 1.d0
   endif

!               Coefficients aux interfaces entre couches
   cnurho = cnu(1)*rho(1)
   cgarho = cgam(1)*rho(1)
   cnugam = cnu(1)*cgam(1)
   cb2 = kr2/ckb2(1)
   ckb2i2 = 1.d0/ckb2(1)

   do ic = 2, nc

      ic1 = ic - 1
      ckb2i1 = ckb2i2
      ckb2i2 = 1.d0/ckb2(ic)
      cb1 = cb2
      cb2 = kr2*ckb2i2
      ca1d = rho(ic1)*(1.d0-2.d0*cb1)
      ca2d = rho(ic)*(1.d0-2.d0*cb2)
      ca = ca2d-ca1d
      cb = ca2d+2.*rho(ic1)*cb1
      cc = ca1d+2.*rho(ic)*cb2
      cd = 2.*(rho(ic)*ckb2i2 - rho(ic1)*ckb2i1)
      ce = cb*cnu(ic1) + cc*cnu(ic)
      cf = cb*cgam(ic1) + cc*cgam(ic)
      cg = ca - cd*cnu(ic1)*cgam(ic)
      cgkr2 = cg*kr2
      ch = ca - cd*cnu(ic)*cgam(ic1)
      chkr2 = ch*kr2
      cdd = 1.d0/(ce*cf + cg*chkr2)
      cdd2 = 2.d0*cdd

      ctmp2 = cnurho*cdd2
      ctmp3 = cgarho*cdd2
      ctmp4 = (ca*cc + cb*cd*cnugam)*cdd2
      cnurho = cnu(ic)*rho(ic)
      cgarho = cgam(ic)*rho(ic)
      cnugam = cnu(ic)*cgam(ic)
      ctmp1 = (ca*cb + cc*cd*cnugam)*cdd2
      ctmp5 = cnurho*cdd2
      ctmp6 = cgarho*cdd2
      ctmp7 = cf*(cb*cnu(ic1) - cc*cnu(ic))
      ctmp8 = chkr2*(ca + cd*cnu(ic1)*cgam(ic))
      ctmp9 = ce*(cb*cgam(ic1) - cc*cgam(ic))
      ctmp10 = cgkr2*(ca + cd*cnu(ic)*cgam(ic1))

      rd(ic, 1, 1) = (ctmp7 - ctmp8)*cdd
      rd(ic, 1, 2) = -kr2*cgam(ic1)*ctmp1*ai*aki
      rd(ic, 2, 2) = -(ctmp9 - ctmp10)*cdd*aki
      rd(ic, 2, 1) = -cnu(ic1)*ctmp1*coef
      td(ic, 1, 1) = ctmp2*cf
      td(ic, 1, 2) = -ctmp3*cgkr2*ai*aki
      td(ic, 2, 2) = ctmp3*ce
      td(ic, 2, 1) = ctmp2*ch*coef*aki

      ru(ic, 1, 1) = -(ctmp7 + ctmp10)*cdd
      ru(ic, 1, 2) = kr2*cgam(ic)*ctmp4*ai
      ru(ic, 2, 2) = (ctmp9 + ctmp8)*cdd*aki
      ru(ic, 2, 1) = cnu(ic)*ctmp4*coef*aki
      tu(ic, 1, 1) = ctmp5*cf
      tu(ic, 1, 2) = ctmp6*chkr2*ai
      tu(ic, 2, 2) = ctmp6*ce
      tu(ic, 2, 1) = -ctmp5*cg*coef

!   Modification pour calculateur a faible dynamique [1.e-300; 1.e+300]
      arg   = -ai*cnu(ic1)*hc(ic1)
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


      cs1 = rho(ic1)*ckb2i1*cgam(ic1)
      cs2 = rho(ic)*ckb2i2*cgam(ic)
      cdelt = 1.d0/(cs1 + cs2)

      rush(ic) = (cs2 - cs1)*cdelt
      rdsh(ic) = -rush(ic)
      tush(ic) = 2.d0*cs2*cdelt
      tdsh(ic) = 2.d0*cs1*cdelt

   end do
!write (6,*) ru(1,:,:),ru(2,:,:)
!write(6,*) rd(1,:,:), rd(2,:,:)
!write(6,*) tu(1,:,:), tu(2,:,:)
!write(6,*) td(1,:,:), td(2,:,:)
!write(6,*)
   return
end
end module
