c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine estrtor3  --  stretch-torsion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "estrtor3" calculates the stretch-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
      subroutine estrtor3
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bndstr
      use bound
      use energi
      use group
      use inform
      use iounit
      use math
      use strtor
      use torpot
      use tors
      use usage
      implicit none
      integer i,k,istrtor
      integer ia,ib,ic,id
      real*8 e,dr
      real*8 fgrp,angle
      real*8 rt2,ru2,rtru
      real*8 rba,rcb,rdc
      real*8 e1,e2,e3
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      logical proceed
      logical header,huge
c
c
c     zero out the stretch-torsion energy and partitioning terms
c
      nebt = 0
      ebt = 0.0d0
      do i = 1, n
         aebt(i) = 0.0d0
      end do
      if (nstrtor .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nstrtor.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Stretch-Torsion Interactions :',
     &           //,' Type',25x,'Atom Names',21x,'Angle',6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nstrtor,ist,itors,kst,bl,
!$OMP& tors1,tors2,tors3,use,x,y,z,storunit,use_group,use_polymer,
!$OMP& name,verbose,debug,header,iout)
!$OMP& shared(ebt,nebt,aebt)
!$OMP DO reduction(+:ebt,nebt,aebt) schedule(guided)
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtor
         i = ist(1,istrtor)
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
            rba = sqrt(xba*xba + yba*yba + zba*zba)
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
            if (min(rba,rcb,rdc) .ne. 0.0d0) then
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
               rt2 = max(rt2,0.000001d0)
               ru2 = xu*xu + yu*yu + zu*zu
               ru2 = max(ru2,0.000001d0)
               rtru = sqrt(rt2*ru2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     compute multiple angle trigonometry and phase terms
c
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     get the stretch-torsion values for the first bond
c
               v1 = kst(1,istrtor)
               v2 = kst(2,istrtor)
               v3 = kst(3,istrtor)
               k = ist(2,istrtor)
               dr = rba - bl(k)
               e1 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the second bond
c
               v1 = kst(4,istrtor)
               v2 = kst(5,istrtor)
               v3 = kst(6,istrtor)
               k = ist(3,istrtor)
               dr = rcb - bl(k)
               e2 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     get the stretch-torsion values for the third bond
c
               v1 = kst(7,istrtor)
               v2 = kst(8,istrtor)
               v3 = kst(9,istrtor)
               k = ist(4,istrtor)
               dr = rdc - bl(k)
               e3 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e1 = e1 * fgrp
                  e2 = e2 * fgrp
                  e3 = e3 * fgrp
               end if

c     increment the total stretch-torsion energy
c
               nebt = nebt + 1
               e = e1 + e2 + e3
               ebt = ebt + e
               aebt(ia) = aebt(ia) + 0.5d0*e1
               aebt(ib) = aebt(ib) + 0.5d0*(e1+e2)
               aebt(ic) = aebt(ic) + 0.5d0*(e2+e3)
               aebt(id) = aebt(id) + 0.5d0*e3
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 3.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Stretch-Torsion',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,30)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   30             format (' StrTors',3x,4(i7,'-',a3),f11.4,f12.4)
               end if
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
