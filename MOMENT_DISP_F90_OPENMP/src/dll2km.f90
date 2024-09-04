!	
!       input x1,x2: array of latitude
!             y1,y2: array of longitude
!       replaced by cartesian coordinates in METERS
!
        subroutine ll2km(x1,y1,n1,x2,y2,n2)
        implicit none
        real*8 x1(*),y1(*),x2(*),y2(*),xmin,xmax,x,y
        integer n1,n2,i

        xmax=x1(1)
        xmin=xmax
        do i=1,n1
        if (x1(i).lt.xmin) then
           xmin=x1(i)
        endif
        if (x1(i).gt.xmax) then
           xmax=x1(i)
        endif
        enddo 
        do i=1,n2
        if (x2(i).lt.xmin) then
           xmin=x2(i)
        endif
        if (x2(i).gt.xmax) then
           xmax=x2(i)
        endif
        enddo 

 	if (xmin == xmax) xmax = xmax+0.01
        call init_lcc(x1(1),y1(1),xmin,xmax)
        do i=1,n1
		call projlcc(x1(i),y1(i),x,y,1)
                x1(i)=y*1000.
                y1(i)=x*1000.
	enddo
        do i=1,n2
		call projlcc(x2(i),y2(i),x,y,1)
                x2(i)=y*1000.
                y2(i)=x*1000.
	enddo

        end

!  SUBROUTINE TO DO LAMBERT CONFORMAL CONIC PROJECTION
!  From USGS Bulletin no.1532, 1982, p. 105
      subroutine projlcc(rlat,rlon,xloc,yloc,isens)

      implicit real*8 (a-h,o-z)
      common /laccproj/ rlon0,radius,p0,fcon,cone,pi4
      common /trig/ pi,radpdeg

!  CONSTANTS SET IN EARLIER ROUTINE
!           cone=alog(cos(radpdeg*stdlat1)/cos(radpdeg*stdlat2))/
!    +         alog(tan(pi4+radpdeg*stdlat2/2.)/
!    +         tan(pi4+radpdeg*stdlat1/2.))
!           fcon=cos(radpdeg*stdlat1)*
!    +         (tan(pi4+radpdeg*stdlat1/2.)**cone)/cone
!           p0=radius*fcon/(tan(pi4+radpdeg*rlat0/2.)**cone)

	if (isens.eq.1) then
      p=radius*fcon/(dtan(pi4+radpdeg*rlat/2.)**cone)
      theta=cone*radpdeg*(rlon-rlon0)
      xloc=p*dsin(theta)
      yloc=p0-p*dcos(theta)
	else
      theta=datan2(xloc, (p0-yloc))
      p=dsqrt(xloc**2 + (p0-yloc)**2)
      rlon = theta/cone/radpdeg+rlon0
      rlat = (datan ((radius*fcon/p)**(1./cone))-pi4)*2./radpdeg
        endif

      return
      end

!  SUBROUTINE TO INITIALIZE LAMBERT CONF CONIC PROJECTION CONSTANTS
      subroutine init_lcc(rlat0,lon0,stdlat1,stdlat2)

      implicit real*8 (a-h,o-z)
      real*8 lon0
      common /laccproj/ rlon0,radius,p0,fcon,cone,pi4
      common /trig/ pi,radpdeg
      rlon0=lon0
      pi4=atan(1.)
      pi=pi4*4.
      radpdeg=pi/180.
      radius=6378.

      cone=dlog(dcos(radpdeg*stdlat1)/dcos(radpdeg*stdlat2))/ &
         dlog(dtan(pi4+radpdeg*stdlat2/2.)/ &
         dtan(pi4+radpdeg*stdlat1/2.)) 
      fcon=dcos(radpdeg*stdlat1)* &
         (dtan(pi4+radpdeg*stdlat1/2.)**cone)/cone
      p0=radius*fcon/(dtan(pi4+radpdeg*rlat0/2.)**cone)

      return
      end
