!**************************************************************************
!                     AXITRA Moment Version
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
!*************************************************************************
module reflect4m
contains

subroutine reflect4(jf, ik, tmin, tconv, nc, nr, ns, ncs, ncr, uflow)
   use dimension1
   use dimension2
   use parameter
   implicit none

   integer :: nc, nr, ns, ncs,ncr
   integer :: jf, ik
   logical :: uflow
! Local
   integer :: ic, ir, ir1, ir2, ir3, idel, is, is1, is2, is3, jrs, ics, is0, ir0
   logical :: tmin, tconv(nr, ns)

   complex(kind=fd) :: egam, enu, s1phiu, s1phid
   complex(kind=fd) :: s1psiu, s1psid, s2phiu, s2phid, s2psiu, s2psid, s3phiu
   complex(kind=fd) :: s3phid, s3psiu, s3psid, s4phid, s4phiu, s4psid, s4psiu
   complex(kind=fd) :: s5, s6, arg, enuinv, egaminv, cs2, cs4, cs3, cs5, cs6, cs7, cs8, cs9
   complex(kind=fd) :: cu, cdu1, cdu2, cu2, cdu3, cdu4, cdu5, cu3, cdu6, cdu7, cdu8, cdu9, cdu10, cdu11
   complex(kind=fd) :: cdu(nr, ns, 11)
   complex(kind=fd) :: cz1, cz1b, cz2, cz2b, cz3, cz3b, cz4, cz4b, cr1, cr2, cr3
   real(kind=fd)     :: zc, dz, r1, i1, r2, i2, zsc

   complex(kind=fd) :: pu(2, 4), pd(2, 4), push(2), pdsh(2)
   complex(kind=fd) :: ftdosh(nc), ftup(nc, 2, 2), ftdo(nc, 2, 2), ftupsh(nc)
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
!   pu(n) = [fup(n)]*[fup(n+1)]* *[fup(ics-1] . su
!
!   d ou l on tire pd(n) par  pd(n) = [nt(n)] . pu(n)
!
!  Couche (m) au dessous de la couche source :
!
!   pd(m) = [fdo(m)]*[fdo(m-1)]* *[fdo(ics+1)] . sd
!
!   d ou l on tire pu(m) par  pu(m) = [mt(m)] . pd(m)
!
!                -------------------------
!   On pose :
!
!        [ftup(n)] = [fup(n)]*...*[fup(ics-1)]*[tud]
!
!        [ftdo(m)] = [fdo(m)]*...*[fdo(ics+1)]*[tdu]
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do is1 = 1, ncs ! boucle sur les couches sources
      ics = isc(is1)

!                Couches au dessus de la couche source
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
!               Constantes dependant du nombre d onde et
!               de la couche source

      cs2 = c2(ics)
      cs4 = kr2
      cs5 = ai*cgam(ics)
      cs3 = ckb2(ics)/cs5
      cs6 = ckb2(ics)
      cs8 = ai*cnu(ics)
      cs7 = cs8*kr2
      cs9 = (ckb2(ics) - 2.*kr2)/cs5

      do ir1 = 1, ncr !boucle sur les couches recepteur
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
!  Deplacements pour chaque sources du tenseur, exprime a l aide de
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

!                          Termes : C5(kr,z)
!                          -----------------

!                Dephasage par rapport au sommet de la couche pour
!                les ondes PHI, PSI et KHI

               arg = -ai*cgam(ic)*zc
               if (dreal(arg) .lt. explim)  then
                   uflow=.true.
                   arg=cmplx(explim,dimag(arg)) 
               endif
               egam = exp(arg)
               egaminv = 1./egam
               arg = -ai*cnu(ic)*zc
               if (dreal(arg) .lt. explim)  then
                   uflow=.true.
                   arg=cmplx(explim,dimag(arg)) 
               endif
               enu=exp(arg)
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

!                        Termes u=C2(kr)*C4(kr.r)*C5(kr,z)
!                       -------------------------------
!                    (Tous les termes de deplacement presentant les memes
!                      dependances en theta sont regroupes)

!                Mxx :
!               -----
!                        PHI0 + PSI0
!        ur => u(1) ; u(2)
!        ut => -u(2)
!        uz => u(3) ; u(4)
!                        KHI0
!        ur => u(5)
!        ut => -u(6)
                        cu = cs2*cz1 + cz4
                        cdu1 = k3(jrs)*cu
                        cdu2 = k4(jrs)*cu
                        u(ir, is, 1) = cdu1 + u(ir, is, 1)
                        u(ir, is, 2) = cdu2 + u(ir, is, 2)
                        cu2 = cs2*cz1b+cz4b
                        cdu3 = k1(jrs)*cu2
                        cdu4 = k2(jrs)*cu2
                        cdu5 = cs3*k4(jrs)*s5
                        cdu6 = cs3*k3(jrs)*s5
                        u(ir, is, 3) = cdu3 + u(ir, is, 3)
                        u(ir, is, 4) = cdu4 + u(ir, is, 4)
                        u(ir, is, 5) = cdu5 + u(ir, is, 5)
                        u(ir, is, 6) = cdu6 + u(ir, is, 6)

!        Double couple Mxy+Myx :
!       -----------------------
!                        PHI0 + PSI0
!             ur => u(1)
!             ut => u(2)
!             uz => u(3)
!                        KHI0
!             ur => u(5)
!             ut => u(6)

!        Double couple Mxz+Mzx :
!       -----------------------
!                        PHI0 + PSI0 + KHI0
!        ur =>  u(7)
!        ut => -u(8)
!        uz =>  u(9)

                        cu3 = -2.*cs4*cz2 + cs9*cz3
                        cdu7 = k5(jrs)*cu3 - cs6*k2(jrs)*s6
                        cdu8 = k2(jrs)*cu3 - cs6*k5(jrs)*s6
                        cdu9 = fj1(jrs)*(-2.*cs4*cz2b+cs9*cz3b)
                        u(ir, is, 7) = cdu7 + u(ir, is, 7)
                        u(ir, is, 8) = cdu8 + u(ir, is, 8)
                        u(ir, is, 9) = cdu9 + u(ir, is, 9)

!                Myy :
!               -----
!                        PHI0 + PSI0
!     ur => u(1) ; u(2)
!     ut => u(2)
!     uz => u(3) ; u(4)
!                        KHI0
!     ur => -u(5)
!     ut => u(6)

!        Double couple Myz+Mzy :
!       -----------------------
!                        PHI0 + PSI0 + KHI0
!     ur =>  u(7)
!     ut =>  u(8)
!     uz =>  u(9)

!                Mzz :
!               -----
!                        PHI0 + PSI0
!        ur => u(10)
!        uz => -u(11)

                        cdu10 = fj1(jrs)*(cs7*cz1 + kr2*cz4)
                        cdu11 = k0(jrs)*(cs8*cz1b+cz4b)
                        u(ir, is, 10) = cdu10+u(ir, is, 10)
                        u(ir, is, 11) = cdu11+u(ir, is, 11)

                        cdu(ir, is, 1) = cdu1
                        cdu(ir, is, 2) = cdu2
                        cdu(ir, is, 3) = cdu3
                        cdu(ir, is, 4) = cdu4
                        cdu(ir, is, 5) = cdu5
                        cdu(ir, is, 6) = cdu6
                        cdu(ir, is, 7) = cdu7
                        cdu(ir, is, 8) = cdu8
                        cdu(ir, is, 9) = cdu9
                        cdu(ir, is, 10) = cdu10
                        cdu(ir, is, 11) = cdu11

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
               cdu6 = cdu(ir, is, 6)
               cdu7 = cdu(ir, is, 7)
               cdu8 = cdu(ir, is, 8)
               cdu9 = cdu(ir, is, 9)
               cdu10 = cdu(ir, is, 10)
               cdu11 = cdu(ir, is, 11)

               r1 = dble(u(ir, is, 1))
               i1 = dimag(u(ir, is, 1))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu1)
               i2 = dimag(cdu1)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 2))
               i1 = dimag(u(ir, is, 2))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu2)
               i2 = dimag(cdu2)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 3))
               i1 = dimag(u(ir, is, 3))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu3)
               i2 = dimag(cdu3)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 4))
               i1 = dimag(u(ir, is, 4))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu4)
               i2 = dimag(cdu4)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 5))
               i1 = dimag(u(ir, is, 5))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu5)
               i2 = dimag(cdu5)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 6))
               i1 = dimag(u(ir, is, 6))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu6)
               i2 = dimag(cdu6)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 7))
               i1 = dimag(u(ir, is, 7))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu7)
               i2 = dimag(cdu7)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 8))
               i1 = dimag(u(ir, is, 8))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu8)
               i2 = dimag(cdu8)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 9))
               i1 = dimag(u(ir, is, 9))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu9)
               i2 = dimag(cdu9)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 10))
               i1 = dimag(u(ir, is, 10))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu10)
               i2 = dimag(cdu10)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               r1 = dble(u(ir, is, 11))
               i1 = dimag(u(ir, is, 11))
               r1 = (r1*r1 + i1*i1)*uconv
               r2 = dble(cdu11)
               i2 = dimag(cdu11)
               r2 = (r2*r2 + i2*i2)
               tconv(ir, is) = ((r2 .le. r1) .and. tconv(ir, is))

               ttconv = ttconv .and. tconv(ir, is)
            endif
         enddo !end loop receiver
      enddo !end loop source
   endif

   return
end
end module
