!******************************************************************************
!                     AXITRA Moment Version
!
!*                     SUBROUTINE INITDATA
!*
!*    Initialisation de divers parametres
!*
!*    Input:
!*        hc,zr,zs,nc,nr
!*    Output:
!*             ncr,irc,nzr,irzz,nzrr,rr
!*             ncs,isc,nzs,iszz,nzss,rr,iss
!*    Modified:
!*           hc
!******************************************************************************
module initdatam
use parameter 
contains

function mean(tab)
implicit none
real(kind=fd) :: mean,tab(:)

integer ::i
  mean=0.d0
  do i=1,size(tab)
     mean=mean+tab(i)
  enddo
  mean=mean/size(tab)
end function mean

subroutine initdata(latlon, nr, ns, nc, ncr, ncs, nrs, rmax)

   use dimension1
   use dimension2
   use parameter
   implicit none
! argument
   logical      :: latlon
! Local
   integer      :: ir, ir1, ir2, ic, jr, jrr, js, jss, is, is1, is2, i, nr, ns, nc
   integer      :: rindex(nr), index(ns), ncs, ncr, nrs
   logical      :: tc
   real(kind=fd) :: hh, tmp, r(nr, ns), rmax
   real(kind=fd) :: zs2(ns)

!++++++++++++
!        read stations and source coordinates
!++++++++++++
   do is = 1, ns
      read (in2, *) index(is), xs(is), ys(is), zs(is)
   enddo
   do ir = 1, nr
      read (in3, *) rindex(ir), xr(ir), yr(ir), zr(ir)
   enddo


!++++++++++++
!        conversion interface -> epaisseur des couches
!++++++++++++

   if (doubleEquality(hc(1),0._8)) then
      do ic = 1, nc - 1
         hc(ic) = hc(ic + 1) - hc(ic)
      enddo
   endif



!++++++++++++
!        on reordonne les sources par profondeur croissante
!++++++++++++

   do is1 = 1, ns - 1
      do is2 = is1, ns
         if (zs(is1) .gt. zs(is2)) then
            tmp = xs(is1)
            xs(is1) = xs(is2)
            xs(is2) = tmp
            tmp = ys(is1)
            ys(is1) = ys(is2)
            ys(is2) = tmp
            tmp = zs(is1)
            zs(is1) = zs(is2)
            zs(is2) = tmp
            tmp = index(is1)
            index(is1) = index(is2)
            index(is2) = tmp
         endif
      enddo
   enddo
   close (in2)
   zs2=zs


!++++++++++++
!       on calcule :
!       ncs: nombre de couches contenant un source
!       isc(): liste des couches contenant un source
!       nzs(i): nbre de sources de prof. differente dans couche i
!       nzss(j,i): nbre de sources a la prof j, dans la couche i
!       izss(,j,i): indice dans xs(),ys(),zs() des sources a la prof j
!                   dans la couche i
!++++++++++++
   nzss=0
   do is = 1, ns
!                       compute ic,zc
      ic = 1
      hh = hc(1)
      do while ((zs(is) .gt. hh) .and. (ic .lt. nc))
         zs(is) = zs(is) - hh
         ic = ic + 1
         hh = hc(ic)
      enddo
      cff(is) = 1./rho(ic)
!                       compute isc(),ncs,js
      if (is .eq. 1) then
         isc(1) = ic
         ncs = 1
         js = 1
      else
         is1 = 1
         tc = .true.
         do while (is1 .le. ncs)
            if (ic .eq. isc(is1)) then
               js = is1
               tc = .false.
            endif
            is1 = is1 + 1
         enddo
         if (tc) then
            ncs = ncs + 1
            isc(ncs) = ic
            js = ncs
            nzs(js) = 0
         endif
      endif
!                       compute nzs(),jss
      if (is .eq. 1) then
         nzs(1) = 1
         jss = 1
         tc = .false.
      else
         is2 = 1
         tc = .true.
         do while (is2 .le. nzs(js))
!            if (zs(is) .eq. zs(izss(1, is2, js))) then
            if (doubleEquality(zs(is),zs(izss(1, is2, js)))) then
               jss = is2
               tc = .false.
            endif
            is2 = is2 + 1
         enddo
      endif
      if (tc) then
         nzs(js) = nzs(js) + 1
         jss = nzs(js)
      endif
!                       compute nzss(,),izss(,,)
      nzss(jss, js) = nzss(jss, js) + 1
      izss(nzss(jss, js), jss, js) = is
   enddo

!++++++++++++
!         on reordonne les stations par profondeur croissante
!++++++++++++
   do ir1 = 1, nr - 1
      do ir2 = ir1, nr
         if (zr(ir1) .gt. zr(ir2)) then
            tmp = xr(ir1)
            xr(ir1) = xr(ir2)
            xr(ir2) = tmp
            tmp = yr(ir1)
            yr(ir1) = yr(ir2)
            yr(ir2) = tmp
            tmp = zr(ir1)
            zr(ir1) = zr(ir2)
            zr(ir2) = tmp
            tmp = rindex(ir1)
            rindex(ir1) = rindex(ir2)
            rindex(ir2) = tmp
         endif
      enddo
   enddo
   close (in3)

!++++++++++++
!        on calcule :
!        ncr: nombre de couches contenant un recepteur
!        irc(): liste des couches contenant un recept
!       nzr(i): nbre de recept. de prof. differente dans couche i
!        nzrr(j,i): nbre de recept a la prof j, dans la couche i
!        izrr(,j,i): indice dans xr(),yr(),zr() des recept a la prof j
!                    dans la couche i
!++++++++++++
   nzrr=0
   do ir = 1, nr
!                         compute ic,zc
      ic = 1
      hh = hc(1)
      do while ((zr(ir) .gt. hh) .and. (ic .lt. nc))
         zr(ir) = zr(ir) - hh
         ic = ic + 1
         hh = hc(ic)
      enddo
!                        compute irc(),ncr,jr
      if (ir .eq. 1) then
         irc(1) = ic
         ncr = 1
         jr = 1
      else
         ir1 = 1
         tc = .true.
         do while (ir1 .le. ncr)
            if (ic .eq. irc(ir1)) then
               jr = ir1
               tc = .false.
            endif
            ir1 = ir1 + 1
         enddo
         if (tc) then
            ncr = ncr + 1
            irc(ncr) = ic
            jr = ncr
            nzr(jr) = 0
         endif
      endif
!                        compute nzr(),jrr
      if (ir .eq. 1) then
         nzr(1) = 1
         jrr = 1
         tc = .false.
      else
         ir2 = 1
         tc = .true.
         do while (ir2 .le. nzr(jr))
!           if (zr(ir) .eq. zr(izrr(1, ir2, jr))) then
            if (doubleEquality(zr(ir),zr(izrr(1, ir2, jr)))) then
               jrr = ir2
               tc = .false.
            endif
            ir2 = ir2 + 1
         enddo
      endif
      if (tc) then
         nzr(jr) = nzr(jr) + 1
         jrr = nzr(jr)
      endif
!                        compute nzrr(,),izrr(,,)
      nzrr(jrr, jr) = nzrr(jrr, jr) + 1
      izrr(nzrr(jrr, jr), jrr, jr) = ir
   enddo

!++++++++++
! convert from lat-lon to METER
!++++++++++
   if (latlon) then
      call ll2km(xr, yr, nr, xs, ys, ns)
      open (20, file='source.xyz', form='formatted')
      do i = 1, ns
         write (20,"(I10,3F15.3)") index(i), xs(i), ys(i), zs2(i)
      enddo
      close (20)
      open (20, file='station.xyz', form='formatted')
      do i = 1, nr
         write (20,"(I10,3F15.3)") rindex(i), xr(i), yr(i), zr(i)
      enddo
      close (20)
   endif

!++++++++++++
!         distances radiales / source
!         on ne garde que les distances differentes, stockees dans
!         rr(). tableau d indirection irr().
!++++++++++++
   rmax=0.d0
   nrs = 0 !calcule dist. rad.
   do is = 1, ns
   do ir = 1, nr
      nrs = nrs + 1
      r(ir, is) = sqrt((xr(ir) - xs(is))*(xr(ir) - xs(is)) + &
                       (yr(ir) - ys(is))*(yr(ir) - ys(is)))
      rr(nrs) = r(ir, is)
      rmax=max(rmax,rr(nrs))

   enddo
   enddo

   ir1 = 1 !elimine dist. rad. egales
   do while (ir1 .lt. nrs)
      ir2 = ir1 + 1
      do while (ir2 .le. nrs)
!     if (rr(ir1) .eq. rr(ir2)) then
      if (doubleEquality(rr(ir1),rr(ir2))) then
         rr(ir2) = rr(nrs)
         nrs = nrs - 1
      else
         ir2 = ir2 + 1
      endif
      enddo
      ir1 = ir1 + 1
   enddo

! Tableau d indirection
   do is = 1, ns
      do ir = 1, nr
         do ir2 = 1, nrs
!           if (r(ir, is) .eq. rr(ir2)) irs(ir, is) = ir2
            if (doubleEquality(r(ir, is),rr(ir2))) irs(ir, is) = ir2
         enddo
      enddo
   enddo

! coef azimut.
   do is = 1, ns
      do ir = 1, nr
         if (r(ir, is) .ne. 0.) then
            cosr(ir, is) = (xr(ir) - xs(is))/r(ir, is)
            sinr(ir, is) = (yr(ir) - ys(is))/r(ir, is)
         else
            cosr(ir, is) = 1.
            sinr(ir, is) = 0.
         endif
      enddo
   enddo

   do ic = 1, nc
      vs2(ic) = vs(ic)*vs(ic)
      vp2(ic) = vp(ic)*vp(ic)
   end do

   return
end
end module
