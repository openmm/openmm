c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrtor1  --  stretch-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrtor1" calculates the stretch-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine estrtor1
      use sizes
      use atoms
      use bndstr
      use bound
      use deriv
      use energi
      use group
      use strtor
      use torpot
      use tors
      use usage
      use virial
      implicit none
      integer i,k,istrtor
      integer ia,ib,ic,id
      real*8 e,dr,fgrp
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
      real*8 dphi1,dphi2,dphi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 ddr,dedphi
      real*8 ddrdx,ddrdy,ddrdz
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the stretch-torsion energy and first derivatives
c
      ebt = 0.0d0
      do i = 1, n
         debt(1,i) = 0.0d0
         debt(2,i) = 0.0d0
         debt(3,i) = 0.0d0
      end do
      if (nstrtor .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nstrtor,ist,itors,kst,bl,
!$OMP& tors1,tors2,tors3,use,x,y,z,storunit,use_group,use_polymer)
!$OMP& shared(ebt,debt,vir)
!$OMP DO reduction(+:ebt,debt,vir) schedule(guided)
c
c     calculate the stretch-torsion energy and first derivatives
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
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca)
                  call image (xdb,ydb,zdb)
               end if
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
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
               dphi1 = cosine*s1 - sine*c1
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
c
c     get the stretch-torsion values for the first bond
c
               v1 = kst(1,istrtor)
               v2 = kst(2,istrtor)
               v3 = kst(3,istrtor)
               k = ist(2,istrtor)
               dr = rba - bl(k)
               e1 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddr = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rba
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e1 = e1 * fgrp
                  dedphi = dedphi * fgrp
                  ddr = ddr * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xba * ddr
               ddrdy = yba * ddr
               ddrdz = zba * ddr
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     determine chain rule components for the first bond
c
               dedxia = zcb*dedyt - ycb*dedzt - ddrdx
               dedyia = xcb*dedzt - zcb*dedxt - ddrdy
               dedzia = ycb*dedxt - xcb*dedyt - ddrdz
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu + ddrdx
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu + ddrdy
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu + ddrdz
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     get the stretch-torsion values for the second bond
c
               v1 = kst(4,istrtor)
               v2 = kst(5,istrtor)
               v3 = kst(6,istrtor)
               k = ist(3,istrtor)
               dr = rcb - bl(k)
               e2 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddr = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rcb
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e2 = e2 * fgrp
                  dedphi = dedphi * fgrp
                  ddr = ddr * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     increment chain rule components for the second bond
c
               dedxia = dedxia + zcb*dedyt - ycb*dedzt
               dedyia = dedyia + xcb*dedzt - zcb*dedxt
               dedzia = dedzia + ycb*dedxt - xcb*dedyt
               dedxib = dedxib + yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu - ddrdx
               dedyib = dedyib + zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu - ddrdy
               dedzib = dedzib + xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu - ddrdz
               dedxic = dedxic + zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu + ddrdx
               dedyic = dedyic + xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu + ddrdy
               dedzic = dedzic + yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu + ddrdz
               dedxid = dedxid + zcb*dedyu - ycb*dedzu
               dedyid = dedyid + xcb*dedzu - zcb*dedxu
               dedzid = dedzid + ycb*dedxu - xcb*dedyu
c
c     get the stretch-torsion values for the third bond
c
               v1 = kst(7,istrtor)
               v2 = kst(8,istrtor)
               v3 = kst(9,istrtor)
               k = ist(4,istrtor)
               dr = rdc - bl(k)
               e3 = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddr = storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rdc
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e3 = e3 * fgrp
                  dedphi = dedphi * fgrp
                  ddr = ddr * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xdc * ddr
               ddrdy = ydc * ddr
               ddrdz = zdc * ddr
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     increment chain rule components for the third bond
c
               dedxia = dedxia + zcb*dedyt - ycb*dedzt
               dedyia = dedyia + xcb*dedzt - zcb*dedxt
               dedzia = dedzia + ycb*dedxt - xcb*dedyt
               dedxib = dedxib + yca*dedzt - zca*dedyt + zdc*dedyu
     &                     - ydc*dedzu
               dedyib = dedyib + zca*dedxt - xca*dedzt + xdc*dedzu
     &                     - zdc*dedxu
               dedzib = dedzib + xca*dedyt - yca*dedxt + ydc*dedxu
     &                     - xdc*dedyu
               dedxic = dedxic + zba*dedyt - yba*dedzt + ydb*dedzu
     &                     - zdb*dedyu - ddrdx
               dedyic = dedyic + xba*dedzt - zba*dedxt + zdb*dedxu
     &                     - xdb*dedzu - ddrdy
               dedzic = dedzic + yba*dedxt - xba*dedyt + xdb*dedyu
     &                     - ydb*dedxu - ddrdz
               dedxid = dedxid + zcb*dedyu - ycb*dedzu + ddrdx
               dedyid = dedyid + xcb*dedzu - zcb*dedxu + ddrdy
               dedzid = dedzid + ycb*dedxu - xcb*dedyu + ddrdz
c
c     increment the stretch-torsion energy and gradient
c
               e = e1 + e2 + e3
               ebt = ebt + e
               debt(1,ia) = debt(1,ia) + dedxia
               debt(2,ia) = debt(2,ia) + dedyia
               debt(3,ia) = debt(3,ia) + dedzia
               debt(1,ib) = debt(1,ib) + dedxib
               debt(2,ib) = debt(2,ib) + dedyib
               debt(3,ib) = debt(3,ib) + dedzib
               debt(1,ic) = debt(1,ic) + dedxic
               debt(2,ic) = debt(2,ic) + dedyic
               debt(3,ic) = debt(3,ic) + dedzic
               debt(1,id) = debt(1,id) + dedxid
               debt(2,id) = debt(2,id) + dedyid
               debt(3,id) = debt(3,id) + dedzid
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
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
