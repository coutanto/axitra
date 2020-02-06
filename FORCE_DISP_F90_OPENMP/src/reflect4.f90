!**************************************************************************
!*                     AXITRA Force Version
!*
!*                     SUBROUTINE REFLECT4
!*
!*  Calcul des potentiels et des deplacements a chaque recepteur.
!*  - Matrice de passage des potentiels de la couche source aux couches
!*    recepteur (FTUP et FTDO)
!*  - Calcul des potentiels dans toutes les couches (PU et PD)
!*  - Calcul de 11 deplacements (termes intermediaires) a chaque recepteur
!*    (U)
!*
!
! input:
!	jf= frequency index
!	ik=	wavenumber index
!	tmin= flag to skip or apply convergence criteria
!	tconv= logical convergence criteria
!	nc= number of layer
!	nr= number of receiver
!	ns= nimber of source
!	ncs= number of source layer
!	ncr= number of receiver layer
!**************************************************************************
module reflect4m

contains
subroutine reflect4(jf, ik, tmin, tconv, nc, nr, ns, ncs, ncr)

   use dimension1
   use dimension2
   use parameter

   implicit none

   integer :: nc,nr,ns,ik,jf,ncr,ncs

   logical :: tmin, tconv(nr, ns)
   integer :: is1, ic, is2, is3, ics, is, ir, jrs, ir3,ir0,idel,is0,ir1,ir2
   real(kind=8) :: r1,r2,i1,i2,zc,dz,zsc
   complex(kind=8) :: egam, enu, s1phiu, s1phid, egaminv, enuinv,arg
   complex(kind=8) :: pu(2, 4), pd(2, 4), push(2), pdsh(2)
   complex(kind=8) :: ftdosh(nc), ftup(nc, 2, 2), ftdo(nc, 2, 2), ftupsh(nc)
   complex(kind=8) :: s1psiu, s1psid, s2phiu, s2phid, s2psiu, s2psid, s3phiu
   complex(kind=8) :: s3phid, s3psiu, s3psid, s4phid, s4phiu, s4psid, s4psiu
   complex(kind=8) :: s5, s6, cs2, cs3, cs4, cs5, cs6, cs7, cs8, cs9, cz1, &
                       cz1b, cz2, cz2b, cz3, cz4, cz4b, cz3b
   complex(kind=8) :: cu, cu2, cdu1, cdu2, cdu3, cdu4, cdu5,cdu(nr, ns, 5),cr3,cr1,cr2


!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Matrices de passage des vecteurs potentiels dans la
!  couche source aux vecteurs potentiel dans chaque couche
!
!                    [ftup] et [ftdo]
!
!               ------------------------
!
!  Les vecteurs potentiels pu() et pd() sont obtenus a
!  partir des vecteurs potentiels su() et sd() dans la
!  couche source par :
!
!  Couche (n) au dessus de la couche source :
!
!   pu(n) = [fup(n)]*[fup(n+1)]* *[fup(isc-1] . su
!
!   d'ou l'on tire pd(n) par  pd(n) = [nt(n)] . pu(n)
!
!  Couche (m) au dessous de la couche source :
!
!   pd(m) = [fdo(m)]*[fdo(m-1)]* *[fdo(isc+1)] . sd
!
!   d'ou l'on tire pu(m) par  pu(m) = [mt(m)] . pd(m)
!
!                -------------------------
!   On pose :
!
!        [ftup(n)] = [fup(n)]*...*[fup(isc-1)]*[tud]
!
!        [ftdo(m)] = [fdo(m)]*...*[fdo(isc+1)]*[tdu]
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do is1 = 1, ncs ! boucle sur les couches sources
      ics = isc(is1)
!                Couches au dessus de la couche source
!
      ftup(ics, 1, 1) = 1.
      ftup(ics, 1, 2) = 0.
      ftup(ics, 2, 1) = 0.
      ftup(ics, 2, 2) = 1.
      ftupsh(ics) = 1.

      do ic = ics - 1, 1, -1

         ftup(ic, 1, 1) = fup(ic, 1, 1)*ftup(ic + 1, 1, 1) + fup(ic, 1, 2)*ftup(ic + 1, 2, 1)
         ftup(ic, 1, 2) = fup(ic, 1, 1)*ftup(ic + 1, 1, 2) + fup(ic, 1, 2)*ftup(ic + 1, 2, 2)
         ftup(ic, 2, 1) = fup(ic, 2, 1)*ftup(ic + 1, 1, 1) + fup(ic, 2, 2)*ftup(ic + 1, 2, 1)
         ftup(ic, 2, 2) = fup(ic, 2, 1)*ftup(ic + 1, 1, 2) + fup(ic, 2, 2)*ftup(ic + 1, 2, 2)
         ftupsh(ic) = fupsh(ic)*ftupsh(ic + 1)

      enddo

!                Couches au dessous de la couche source
!
      ftdo(ics, 1, 1) = 1.
      ftdo(ics, 1, 2) = 0.
      ftdo(ics, 2, 1) = 0.
      ftdo(ics, 2, 2) = 1.
      ftdosh(ics) = 1.

      do ic = ics + 1, nc

         ftdo(ic, 1, 1) = fdo(ic, 1, 1)*ftdo(ic - 1, 1, 1) + fdo(ic, 1, 2)*ftdo(ic - 1, 2, 1)
         ftdo(ic, 1, 2) = fdo(ic, 1, 1)*ftdo(ic - 1, 1, 2) + fdo(ic, 1, 2)*ftdo(ic - 1, 2, 2)
         ftdo(ic, 2, 1) = fdo(ic, 2, 1)*ftdo(ic - 1, 1, 1) + fdo(ic, 2, 2)*ftdo(ic - 1, 2, 1)
         ftdo(ic, 2, 2) = fdo(ic, 2, 1)*ftdo(ic - 1, 1, 2) + fdo(ic, 2, 2)*ftdo(ic - 1, 2, 2)
         ftdosh(ic) = fdosh(ic)*ftdosh(ic - 1)

      enddo

!                          Termes : C2(kr)
!                          ---------------
!                Constantes dependant du nombre d'onde et
!                de la couche source

      cs2 = c2(ics)
      cs4 = kr2
      cs5 = ai*cgam(ics)
      cs3 = ckb2(ics)/cs5
      cs6 = ckb2(ics)
      cs8 = ai*cnu(ics)
      cs7 = cs8*kr2
      cs9 = (ckb2(ics) - 2.*kr2)/cs5

      do ir1 = 1, ncr
         ic = irc(ir1)
         cr1 = ai*cgam(ic)
         cr2 = cgam(ic)*cgam(ic)
         cr3 = ai*cnu(ic)
         idel = ics - ic

         do is2 = 1, nzs(is1) !boucle sur les prof. sources dans ics
            is0 = izss(1, is2, is1)
            zsc = zs(is0)

            do ir2 = 1, nzr(ir1) !boucle sur les prof. recept. dans ic=irc(ir1)
               ir0 = izrr(1, ir2, ir1)
               zc = zr(ir0)

               if (idel .eq. 0) then
                  dz = zsc - zc
               else
                  dz = idel
               endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Vecteurs potentiel montant (pu) et descendant (pd),
!  dans chaque couche recepteur, pour les 6 sources elementaires
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                Recepteurs au dessus de la source
               if (dz .gt. 0.) then
                  pu(1, 1) = ftup(ic, 1, 1)*su1(is0, 1) + ftup(ic, 1, 2)*su1(is0, 2)
                  pu(2, 1) = ftup(ic, 2, 1)*su1(is0, 1) + ftup(ic, 2, 2)*su1(is0, 2)
                  pd(1, 1) = nt(ic, 1, 1)*pu(1, 1) + nt(ic, 1, 2)*pu(2, 1)
                  pd(2, 1) = nt(ic, 2, 1)*pu(1, 1) + nt(ic, 2, 2)*pu(2, 1)

                  pu(1, 2) = ftup(ic, 1, 1)*su2(is0, 1) + ftup(ic, 1, 2)*su2(is0, 2)
                  pu(2, 2) = ftup(ic, 2, 1)*su2(is0, 1) + ftup(ic, 2, 2)*su2(is0, 2)
                  pd(1, 2) = nt(ic, 1, 1)*pu(1, 2) + nt(ic, 1, 2)*pu(2, 2)
                  pd(2, 2) = nt(ic, 2, 1)*pu(1, 2) + nt(ic, 2, 2)*pu(2, 2)

                  pu(1, 3) = ftup(ic, 1, 1)*su3(is0, 1) + ftup(ic, 1, 2)*su3(is0, 2)
                  pu(2, 3) = ftup(ic, 2, 1)*su3(is0, 1) + ftup(ic, 2, 2)*su3(is0, 2)
                  pd(1, 3) = nt(ic, 1, 1)*pu(1, 3) + nt(ic, 1, 2)*pu(2, 3)
                  pd(2, 3) = nt(ic, 2, 1)*pu(1, 3) + nt(ic, 2, 2)*pu(2, 3)

                  pu(1, 4) = ftup(ic, 1, 1)*su4(is0, 1) + ftup(ic, 1, 2)*su4(is0, 2)
                  pu(2, 4) = ftup(ic, 2, 1)*su4(is0, 1) + ftup(ic, 2, 2)*su4(is0, 2)
                  pd(1, 4) = nt(ic, 1, 1)*pu(1, 4) + nt(ic, 1, 2)*pu(2, 4)
                  pd(2, 4) = nt(ic, 2, 1)*pu(1, 4) + nt(ic, 2, 2)*pu(2, 4)

                  push(1) = ftupsh(ic)*su1sh(is0)
                  pdsh(1) = ntsh(ic)*push(1)

                  push(2) = ftupsh(ic)*su2sh(is0)
                  pdsh(2) = ntsh(ic)*push(2)

!             Recepteurs au dessous de la source
               else if (dz .lt. 0.) then
                  pd(1, 1) = ftdo(ic, 1, 1)*sd1(is0, 1) + ftdo(ic, 1, 2)*sd1(is0, 2)
                  pd(2, 1) = ftdo(ic, 2, 1)*sd1(is0, 1) + ftdo(ic, 2, 2)*sd1(is0, 2)
                  pu(1, 1) = mt(ic, 1, 1)*pd(1, 1) + mt(ic, 1, 2)*pd(2, 1)
                  pu(2, 1) = mt(ic, 2, 1)*pd(1, 1) + mt(ic, 2, 2)*pd(2, 1)

                  pd(1, 2) = ftdo(ic, 1, 1)*sd2(is0, 1) + ftdo(ic, 1, 2)*sd2(is0, 2)
                  pd(2, 2) = ftdo(ic, 2, 1)*sd2(is0, 1) + ftdo(ic, 2, 2)*sd2(is0, 2)
                  pu(1, 2) = mt(ic, 1, 1)*pd(1, 2) + mt(ic, 1, 2)*pd(2, 2)
                  pu(2, 2) = mt(ic, 2, 1)*pd(1, 2) + mt(ic, 2, 2)*pd(2, 2)

                  pd(1, 3) = ftdo(ic, 1, 1)*sd3(is0, 1) + ftdo(ic, 1, 2)*sd3(is0, 2)
                  pd(2, 3) = ftdo(ic, 2, 1)*sd3(is0, 1) + ftdo(ic, 2, 2)*sd3(is0, 2)
                  pu(1, 3) = mt(ic, 1, 1)*pd(1, 3) + mt(ic, 1, 2)*pd(2, 3)
                  pu(2, 3) = mt(ic, 2, 1)*pd(1, 3) + mt(ic, 2, 2)*pd(2, 3)

                  pd(1, 4) = ftdo(ic, 1, 1)*sd4(is0, 1) + ftdo(ic, 1, 2)*sd4(is0, 2)
                  pd(2, 4) = ftdo(ic, 2, 1)*sd4(is0, 1) + ftdo(ic, 2, 2)*sd4(is0, 2)
                  pu(1, 4) = mt(ic, 1, 1)*pd(1, 4) + mt(ic, 1, 2)*pd(2, 4)
                  pu(2, 4) = mt(ic, 2, 1)*pd(1, 4) + mt(ic, 2, 2)*pd(2, 4)

                  pdsh(1) = ftdosh(ic)*sd1sh(is0)
                  push(1) = mtsh(ic)*pdsh(1)

                  pdsh(2) = ftdosh(ic)*sd2sh(is0)
                  push(2) = mtsh(ic)*pdsh(2)
               endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Deplacements pour chaque sources du tenseur, exprime a l'aide de
!  de source intermediaires. Chaque source intermediaire correspond
!  aux rayonnements des trois potentiels PHI, PSI, KHI de chaque
!  moment du tenseur.
!
!  ex : Mxy -> PHI0, PSI0, KHI0
!
!              PHI0 -> PHI, PSI apres conversion sur les interfaces
!              PSI0 -> PHI, PSI       "                "
!              KHI0 -> KHI            "                "
!
!                       -------------------------
!
!                u = C2(kr)*C4(kr.r)*C5(kr,z)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!               ic = irc(ir) ??

!                          Termes : C5(kr,z)
!                          -----------------

!                Dephasage par rapport au sommet de la couche pour
!                les ondes PHI, PSI et KHI

               arg = -ai*cgam(ic)*zc
               if (dble(arg) .lt. explim) arg = cmplx(explim, imag(arg))
               egam = exp(arg)
               egaminv = 1./egam
               arg = -ai*cnu(ic)*zc
               if (dble(arg) .lt. explim) arg = cmplx(explim, imag(arg))
               enu = exp(arg)
               enuinv = 1./enu

!                termes sources
               s1phiu = pu(1, 1)*enuinv
               s1phid = pd(1, 1)*enu
               s1psiu = pu(2, 1)*egaminv
               s1psid = pd(2, 1)*egam
               s2phiu = pu(1, 2)*enuinv
               s2phid = pd(1, 2)*enu
               s2psiu = pu(2, 2)*egaminv
               s2psid = pd(2, 2)*egam
               s3phiu = pu(1, 3)*enuinv
               s3phid = pd(1, 3)*enu
               s3psiu = pu(2, 3)*egaminv
               s3psid = pd(2, 3)*egam
               s4phiu = pu(1, 4)*enuinv
               s4phid = pd(1, 4)*enu
               s4psiu = pu(2, 4)*egaminv
               s4psid = pd(2, 4)*egam
               s5 = push(1)*egaminv + pdsh(1)*egam
               s6 = push(2)*egaminv + pdsh(2)*egam

!                                Source phi
               cz1 = (s1phiu + s1phid) + cr1*(s1psiu - s1psid)
               cz1b = cr3*(s1phiu - s1phid) + cs4*(s1psiu + s1psid)
!                                Source phi*sign(z-z0)
               cz2 = (s2phiu + s2phid) + cr1*(s2psiu - s2psid)
               cz2b = cr3*(s2phiu - s2phid) + cs4*(s2psiu + s2psid)
!                                Source psi
               cz3 = (s3phiu + s3phid) + cr1*(s3psiu - s3psid)
               cz3b = cr3*(s3phiu - s3phid) + cs4*(s3psiu + s3psid)
!                                Source psi*sign(z-z0)
               cz4 = (s4phiu + s4phid) + cr1*(s4psiu - s4psid)
               cz4b = cr3*(s4phiu - s4phid) + cs4*(s4psiu + s4psid)

               do is3 = 1, nzss(is2, is1) !boucle distances radiales
                  is = izss(is3, is2, is1)

                  do ir3 = 1, nzrr(ir2, ir1) !boucle sur distances radiales
                     ir = izrr(ir3, ir2, ir1)
                     jrs = irs(ir, is)

                     if (.not. tconv(ir, is)) then !on a pas encore converge
!        Fx et Fy composantes horizontales de la source
!        ur => u(1)
!        ut => u(2)
!        uz => u(3)
#define MV08052018
                        cu = cs2*cz1 + cz4
                        cu2 = cs3*s5
                        cdu1 = k5(jrs)*cu + k2(jrs)*cu2
#ifdef MV08052018
                        cdu2 =  (k2(jrs)*cu + k5(jrs)*cu2)
#else
                        cdu2 = -(k2(jrs)*cu + k5(jrs)*cu2)
#endif
                        cdu3 = j1(jrs)*(cs2*cz1b+cz4b)
                        u(ir, is, 1) = cdu1 + u(ir, is, 1)
                        u(ir, is, 2) = cdu2 + u(ir, is, 2)
                        u(ir, is, 3) = cdu3 + u(ir, is, 3)

!        Fz composante verticale de la source
!        ur => u(4)
!        ut = 0
!        uz => u(5)

                        cdu4 = -j1(jrs)*cs4*(cz2 + cz3/cs5)
                        cdu5 = j0(jrs)*kr*(cz2b+cz3b/cs5)
                        u(ir, is, 4) = cdu4 + u(ir, is, 4)
                        u(ir, is, 5) = cdu5 + u(ir, is, 5)

                        cdu(ir, is, 1) = cdu1
                        cdu(ir, is, 2) = cdu2
                        cdu(ir, is, 3) = cdu3
                        cdu(ir, is, 4) = cdu4
                        cdu(ir, is, 5) = cdu5

                     endif ! test conv deja obtenue
                  enddo ! boucle dist. radiale
               enddo ! boucle dist. radiale
            enddo ! boucle prof. dans couche ic
         enddo ! boucle prof. dans couche ics
      enddo ! boucle sur couche ic
   enddo ! boucle sur couche ics

! Do not try to converge before kmin
   if (tmin) then
      ttconv = .true.
      do is = 1, ns
         do ir = 1, nr
            if (.not. tconv(ir, is)) then
               tconv(ir, is) = .true.
               cdu1 = cdu(ir, is, 1)
               cdu2 = cdu(ir, is, 2)
               cdu3 = cdu(ir, is, 3)
               cdu4 = cdu(ir, is, 4)
               cdu5 = cdu(ir, is, 5)

               r1 = dble(u(ir, is, 1))
               i1 = imag(u(ir, is, 1))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu1)
               i2 = imag(cdu1)
               r2 = (r2*r2 + i2*i2)
!write(6,*) r1,r2
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 2))
               i1 = imag(u(ir, is, 2))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu2)
               i2 = imag(cdu2)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))
!write(6,*) r1,r2
               r1 = dble(u(ir, is, 3))
               i1 = imag(u(ir, is, 3))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu3)
               i2 = imag(cdu3)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))
!write(6,*) r1,r2
               r1 = dble(u(ir, is, 4))
               i1 = imag(u(ir, is, 4))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu4)
               i2 = imag(cdu4)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))
!write(6,*) r1,r2
               r1 = dble(u(ir, is, 5))
               i1 = imag(u(ir, is, 5))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu5)
               i2 = imag(cdu5)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))
!write(6,*) r1,r2
               ttconv = ttconv .and. tconv(ir, is)
!        if (ttconv) write(0,*) 'ir,is,ik',ir,is,ik
            endif
         enddo !end loop receiver
      enddo !end loop source
   endif

   return
end

end module
