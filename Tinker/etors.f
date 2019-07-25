c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine etors  --  torsional angle potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "etors" calculates the torsional potential energy
c
c
      subroutine etors
      use sizes
      use warp
      implicit none
c
c
c     choose standard or potential energy smoothing version
c
      if (use_smooth) then
         call etors0b
      else
         call etors0a
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine etors0a  --  standard torsional angle energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "etors0a" calculates the torsional potential energy
c     using a standard sum of Fourier terms
c
c
      subroutine etors0a
      use sizes
      use atoms
      use bound
      use energi
      use group
      use torpot
      use tors
      use usage
      implicit none
      integer i,ia,ib,ic,id
      real*8 e,rcb,fgrp
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2
      real*8 cosine3,cosine4
      real*8 cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 phi1,phi2,phi3
      real*8 phi4,phi5,phi6
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xdc,ydc,zdc
      real*8 xcb,ycb,zcb
      logical proceed
c
c
c     zero out the torsional potential energy
c
      et = 0.0d0
      if (ntors .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ntors,itors,tors1,tors2,tors3,
!$OMP& tors4,tors5,tors6,use,x,y,z,torsunit,use_group,use_polymer)
!$OMP& shared(et)
!$OMP DO reduction(+:et) schedule(guided)
c
c     calculate the torsional angle energy term
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               phi4 = 1.0d0 + (cosine4*c4 + sine4*s4)
               phi5 = 1.0d0 + (cosine5*c5 + sine5*s5)
               phi6 = 1.0d0 + (cosine6*c6 + sine6*s6)
c
c     calculate the torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                            + v4*phi4 + v5*phi5 + v6*phi6)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               et = et + e
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine etors0b  --  torsional energy for smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "etors0b" calculates the torsional potential energy
c     for use with potential energy smoothing methods
c
c
      subroutine etors0b
      use sizes
      use atoms
      use energi
      use group
      use math
      use torpot
      use tors
      use usage
      use warp
      implicit none
      integer i,ia,ib,ic,id
      real*8 e,rcb,fgrp
      real*8 width,wterm
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2
      real*8 cosine3,cosine4
      real*8 cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 damp1,damp2,damp3
      real*8 damp4,damp5,damp6
      real*8 phi1,phi2,phi3
      real*8 phi4,phi5,phi6
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xdc,ydc,zdc
      real*8 xcb,ycb,zcb
      logical proceed
c
c
c     zero out the torsional potential energy
c
      et = 0.0d0
      if (ntors .eq. 0)  return
c
c     set the extent of smoothing to be performed
c
      width = difft * deform
      if (width .le. 0.0d0) then
         damp1 = 1.0d0
         damp2 = 1.0d0
         damp3 = 1.0d0
         damp4 = 1.0d0
         damp5 = 1.0d0
         damp6 = 1.0d0
      else if (use_dem) then
         damp1 = exp(-width)
         damp2 = exp(-4.0d0*width)
         damp3 = exp(-9.0d0*width)
         damp4 = exp(-16.0d0*width)
         damp5 = exp(-25.0d0*width)
         damp6 = exp(-36.0d0*width)
      else if (use_gda) then
         wterm = difft / 12.0d0
      else if (use_tophat .or. use_stophat) then
         damp1 = 0.0d0
         damp2 = 0.0d0
         damp3 = 0.0d0
         damp4 = 0.0d0
         damp5 = 0.0d0
         damp6 = 0.0d0
         if (width .lt. pi)  damp1 = sin(width) / width
         wterm = 2.0d0 * width
         if (wterm .lt. pi)  damp2 = sin(wterm) / wterm
         wterm = 3.0d0 * width
         if (wterm .lt. pi)  damp3 = sin(wterm) / wterm
         wterm = 4.0d0 * width
         if (wterm .lt. pi)  damp4 = sin(wterm) / wterm
         wterm = 5.0d0 * width
         if (wterm .lt. pi)  damp5 = sin(wterm) / wterm
         wterm = 6.0d0 * width
         if (wterm .lt. pi)  damp6 = sin(wterm) / wterm
      end if
c
c     calculate the torsional angle energy term
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               phi4 = 1.0d0 + (cosine4*c4 + sine4*s4)
               phi5 = 1.0d0 + (cosine5*c5 + sine5*s5)
               phi6 = 1.0d0 + (cosine6*c6 + sine6*s6)
c
c     transform the potential function via smoothing
c
               if (use_gda) then
                  width = wterm * (m2(ia)+m2(ib)+m2(ic)+m2(id))
                  damp1 = exp(-width)
                  damp2 = exp(-4.0d0*width)
                  damp3 = exp(-9.0d0*width)
                  damp4 = exp(-16.0d0*width)
                  damp5 = exp(-25.0d0*width)
                  damp6 = exp(-36.0d0*width)
               end if
               phi1 = phi1 * damp1
               phi2 = phi2 * damp2
               phi3 = phi3 * damp3
               phi4 = phi4 * damp4
               phi5 = phi5 * damp5
               phi6 = phi6 * damp6
c
c     calculate the torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                            + v4*phi4 + v5*phi5 + v6*phi6)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               et = et + e
            end if
         end if
      end do
      return
      end
