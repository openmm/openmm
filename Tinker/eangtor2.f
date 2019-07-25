c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangtor2  --  atomwise angle-torsion Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangtor2" calculates the angle-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine eangtor2 (i)
      use sizes
      use angbnd
      use angtor
      use atoms
      use bound
      use group
      use hessn
      use math
      use torpot
      use tors
      implicit none
      integer i,j,k,iangtor
      integer ia,ib,ic,id
      real*8 dedphi,d2edphi2,fgrp
      real*8 rt2,ru2,rtru,rcb
      real*8 rba2,rcb2,rdc2
      real*8 dot,dt,d2dt
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 terma,termb
      real*8 termc,termd
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 angle,cosang
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xab,yab,zab
      real*8 xbc,ybc,zbc
      real*8 xrab,yrab,zrab
      real*8 xrcb,yrcb,zrcb
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 xrbc,yrbc,zrbc
      real*8 xrdc,yrdc,zrdc
      real*8 xbcp,ybcp,zbcp
      real*8 xdcp,ydcp,zdcp
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 d2phi1,d2phi2,d2phi3
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxia,dyia,dzia
      real*8 dxib,dyib,dzib
      real*8 dxic,dyic,dzic
      real*8 dxid,dyid,dzid
      real*8 dxiaxia,dyiayia,dziazia
      real*8 dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic
      real*8 dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia
      real*8 dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic
      real*8 dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib
      real*8 dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid
      real*8 dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid
      real*8 dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid
      real*8 dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid
      real*8 dzicxid,dzicyid,dziczid
      real*8 domegadxia,domegadyia,domegadzia
      real*8 domegadxib,domegadyib,domegadzib
      real*8 domegadxic,domegadyic,domegadzic
      real*8 domegadxid,domegadyid,domegadzid
      real*8 doxiaxia,doyiayia,doziazia
      real*8 doxibxib,doyibyib,dozibzib
      real*8 doxicxic,doyicyic,doziczic
      real*8 doxidxid,doyidyid,dozidzid
      real*8 doxiayia,doxiazia,doyiazia
      real*8 doxibyib,doxibzib,doyibzib
      real*8 doxicyic,doxiczic,doyiczic
      real*8 doxidyid,doxidzid,doyidzid
      real*8 doxiaxic,doxiayic,doxiazic
      real*8 doyiaxic,doyiayic,doyiazic
      real*8 doziaxic,doziayic,doziazic
      real*8 doxibxic,doxibyic,doxibzic
      real*8 doyibxic,doyibyic,doyibzic
      real*8 dozibxic,dozibyic,dozibzic
      real*8 doxibxid,doxibyid,doxibzid
      real*8 doyibxid,doyibyid,doyibzid
      real*8 dozibxid,dozibyid,dozibzid
      real*8 doxicxid,doxicyid,doxiczid
      real*8 doyicxid,doyicyid,doyiczid
      real*8 dozicxid,dozicyid,doziczid
      real*8 doxibxia,doxibyia,doxibzia
      real*8 doyibxia,doyibyia,doyibzia
      real*8 dozibxia,dozibyia,dozibzia
      real*8 doxicxib,doxicyib,doxiczib
      real*8 doyicxib,doyicyib,doyiczib
      real*8 dozicxib,dozicyib,doziczib
      logical proceed
c
c
c     calculate the angle-torsion interaction Hessian elements
c
      do iangtor = 1, nangtor
         j = iat(1,iangtor)
         ia = itors(1,j)
         ib = itors(2,j)
         ic = itors(3,j)
         id = itors(4,j)
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,0,0)
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
            xab = -xba
            yab = -yba
            zab = -zba
            xbc = -xcb
            ybc = -ycb
            zbc = -zcb
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xab,yab,zab)
               call image (xbc,ybc,zbc)
            end if
            rba2 = xba*xba + yba*yba + zba*zba
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
            if (min(rba2,rcb2,rdc2) .ne. 0.0d0) then
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
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
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
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     set the angle-torsion parameters for the first angle
c
               v1 = kant(1,iangtor)
               v2 = kant(2,iangtor)
               v3 = kant(3,iangtor)
               k = iat(2,iangtor)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosang = dot / sqrt(rba2*rcb2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               dedphi = atorunit * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = atorunit * dt
     &                       * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               d2dt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dt = d2dt * fgrp
               end if
c
c     first and second derivative components for the first angle
c
               terma = -1.0d0 / (rba2*sqrt(rt2))
               termc = 1.0d0 / (rcb2*sqrt(rt2))
               domegadxia = terma * (zba*yt-yba*zt)
               domegadyia = terma * (xba*zt-zba*xt)
               domegadzia = terma * (yba*xt-xba*yt)
               domegadxic = termc * (ycb*zt-zcb*yt)
               domegadyic = termc * (zcb*xt-xcb*zt)
               domegadzic = termc * (xcb*yt-ycb*xt)
               domegadxib = -domegadxia - domegadxic
               domegadyib = -domegadyia - domegadyic
               domegadzib = -domegadzia - domegadzic
c
c     abbreviations used in defining chain rule terms
c
               xrab = 2.0d0 * xab / rba2
               yrab = 2.0d0 * yab / rba2
               zrab = 2.0d0 * zab / rba2
               xrcb = 2.0d0 * xcb / rcb2
               yrcb = 2.0d0 * ycb / rcb2
               zrcb = 2.0d0 * zcb / rcb2
               xabp = (yab*zt-zab*yt) / rt2
               yabp = (zab*xt-xab*zt) / rt2
               zabp = (xab*yt-yab*xt) / rt2
               xcbp = (ycb*zt-zcb*yt) / rt2
               ycbp = (zcb*xt-xcb*zt) / rt2
               zcbp = (xcb*yt-ycb*xt) / rt2
c
c     chain rule terms for second derivative components
c
               doxiaxia = terma*(xab*xcb-dot) + domegadxia*(xcbp-xrab)
               doxiayia = terma*(zt+yab*xcb) + domegadxia*(ycbp-yrab)
               doxiazia = terma*(zab*xcb-yt) + domegadxia*(zcbp-zrab)
               doyiayia = terma*(yab*ycb-dot) + domegadyia*(ycbp-yrab)
               doyiazia = terma*(xt+zab*ycb) + domegadyia*(zcbp-zrab)
               doziazia = terma*(zab*zcb-dot) + domegadzia*(zcbp-zrab)
               doxicxic = termc*(dot-xab*xcb) - domegadxic*(xabp+xrcb)
               doxicyic = termc*(zt-ycb*xab) - domegadxic*(yabp+yrcb)
               doxiczic = -termc*(yt+zcb*xab) - domegadxic*(zabp+zrcb)
               doyicyic = termc*(dot-yab*ycb) - domegadyic*(yabp+yrcb)
               doyiczic = termc*(xt-zcb*yab) - domegadyic*(zabp+zrcb)
               doziczic = termc*(dot-zab*zcb) - domegadzic*(zabp+zrcb)
               doxiaxic = terma*(yab*yab+zab*zab) - domegadxia*xabp
               doxiayic = -terma*xab*yab - domegadxia*yabp
               doxiazic = -terma*xab*zab - domegadxia*zabp
               doyiaxic = -terma*xab*yab - domegadyia*xabp
               doyiayic = terma*(xab*xab+zab*zab) - domegadyia*yabp
               doyiazic = -terma*yab*zab - domegadyia*zabp
               doziaxic = -terma*xab*zab - domegadzia*xabp
               doziayic = -terma*yab*zab - domegadzia*yabp
               doziazic = terma*(xab*xab+yab*yab) - domegadzia*zabp
c
c     get some second derivative chain rule terms by difference
c
               doxibxia = -doxiaxia - doxiaxic
               doxibyia = -doxiayia - doyiaxic
               doxibzia = -doxiazia - doziaxic
               doyibxia = -doxiayia - doxiayic
               doyibyia = -doyiayia - doyiayic
               doyibzia = -doyiazia - doziayic
               dozibxia = -doxiazia - doxiazic
               dozibyia = -doyiazia - doyiazic
               dozibzia = -doziazia - doziazic
               doxibxic = -doxicxic - doxiaxic
               doxibyic = -doxicyic - doxiayic
               doxibzic = -doxiczic - doxiazic
               doyibxic = -doxicyic - doyiaxic
               doyibyic = -doyicyic - doyiayic
               doyibzic = -doyiczic - doyiazic
               dozibxic = -doxiczic - doziaxic
               dozibyic = -doyiczic - doziayic
               dozibzic = -doziczic - doziazic
               doxibxib = -doxibxia - doxibxic
               doxibyib = -doxibyia - doxibyic
               doxibzib = -doxibzia - doxibzic
               doyibyib = -doyibyia - doyibyic
               doyibzib = -doyibzia - doyibzic
               dozibzib = -dozibzia - dozibzic
c
c     scale the first derivatives of the first angle
c
               domegadxia = domegadxia * radian
               domegadyia = domegadyia * radian
               domegadzia = domegadzia * radian
               domegadxic = domegadxic * radian
               domegadyic = domegadyic * radian
               domegadzic = domegadzic * radian
               domegadxib = domegadxib * radian
               domegadyib = domegadyib * radian
               domegadzib = domegadzib * radian
c
c     scale the second derivatives of the first angle
c
               doxiaxia = doxiaxia * d2dt
               doxiayia = doxiayia * d2dt
               doxiazia = doxiazia * d2dt
               doyiayia = doyiayia * d2dt
               doyiazia = doyiazia * d2dt
               doziazia = doziazia * d2dt
               doxicxic = doxicxic * d2dt
               doxicyic = doxicyic * d2dt
               doxiczic = doxiczic * d2dt
               doyicyic = doyicyic * d2dt
               doyiczic = doyiczic * d2dt
               doziczic = doziczic * d2dt
               doxiaxic = doxiaxic * d2dt
               doxiayic = doxiayic * d2dt
               doxiazic = doxiazic * d2dt
               doyiaxic = doyiaxic * d2dt
               doyiayic = doyiayic * d2dt
               doyiazic = doyiazic * d2dt
               doziaxic = doziaxic * d2dt
               doziayic = doziayic * d2dt
               doziazic = doziazic * d2dt
               doxibxia = doxibxia * d2dt
               doxibyia = doxibyia * d2dt
               doxibzia = doxibzia * d2dt
               doyibxia = doyibxia * d2dt
               doyibyia = doyibyia * d2dt
               doyibzia = doyibzia * d2dt
               dozibxia = dozibxia * d2dt
               dozibyia = dozibyia * d2dt
               dozibzia = dozibzia * d2dt
               doxibxic = doxibxic * d2dt
               doxibyic = doxibyic * d2dt
               doxibzic = doxibzic * d2dt
               doyibxic = doyibxic * d2dt
               doyibyic = doyibyic * d2dt
               doyibzic = doyibzic * d2dt
               dozibxic = dozibxic * d2dt
               dozibyic = dozibyic * d2dt
               dozibzic = dozibzic * d2dt
               doxibxib = doxibxib * d2dt
               doxibyib = doxibyib * d2dt
               doxibzib = doxibzib * d2dt
               doyibyib = doyibyib * d2dt
               doyibzib = doyibzib * d2dt
               dozibzib = dozibzib * d2dt
c
c     abbreviations for first derivative chain rule terms
c
               dphidxt = (yt*zcb-ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb-zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb-xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb-ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb-zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb-xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     intermediate terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for first derivative components
c
               dxia = dedphi * dphidxia
               dyia = dedphi * dphidyia
               dzia = dedphi * dphidzia
               dxib = dedphi * dphidxib
               dyib = dedphi * dphidyib
               dzib = dedphi * dphidzib
               dxic = dedphi * dphidxic
               dyic = dedphi * dphidyic
               dzic = dedphi * dphidzic
               dxid = dedphi * dphidxid
               dyid = dedphi * dphidyid
               dzid = dedphi * dphidzid
               dedphi = dedphi * dt
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxiaxib = -dxiaxia - dxiaxic - dxiaxid
               dxiayib = -dxiayia - dxiayic - dxiayid
               dxiazib = -dxiazia - dxiazic - dxiazid
               dyiayib = -dyiayia - dyiayic - dyiayid
               dyiazib = -dyiazia - dyiazic - dyiazid
               dziazib = -dziazia - dziazic - dziazid
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             + 2.0d0 * domegadxia*dxia
     &                             + doxiaxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + domegadxia*dyia + domegadyia*dxia
     &                             + doxiayia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + domegadxia*dzia + domegadzia*dxia
     &                             + doxiazia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + domegadyia*dxia + domegadxia*dyia
     &                             + doxiayia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             + 2.0d0 * domegadyia * dyia
     &                             + doyiayia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + domegadyia*dzia + domegadzia*dyia
     &                             + doyiazia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + domegadxia*dzia + domegadzia*dxia
     &                             + doxiazia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + domegadyia*dzia + domegadzia*dyia
     &                             + doyiazia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             + 2.0d0 * domegadzia * dzia
     &                             + doziazia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + domegadxia*dxib + domegadxib*dxia
     &                             + doxibxia
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + domegadyia*dxib + domegadxib*dyia
     &                             + doxibyia
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + domegadzia*dxib + domegadxib*dzia
     &                             + doxibzia
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + domegadxia*dyib + domegadyib*dxia
     &                             + doyibxia
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + domegadyia*dyib + domegadyib*dyia
     &                             + doyibyia
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + domegadzia*dyib + domegadyib*dzia
     &                             + doyibzia
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + domegadxia*dzib + domegadzib*dxia
     &                             + dozibxia
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + domegadyia*dzib + domegadzib*dyia
     &                             + dozibyia
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + domegadzia*dzib + domegadzib*dzia
     &                             + dozibzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + domegadxia*dxic + domegadxic*dxia
     &                             + doxiaxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + domegadyia*dxic + domegadxic*dyia
     &                             + doyiaxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + domegadzia*dxic + domegadxic*dzia
     &                             + doziaxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + domegadxia*dyic + domegadyic*dxia
     &                             + doxiayic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + domegadyia*dyic + domegadyic*dyia
     &                             + doyiayic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + domegadzia*dyic + domegadyic*dzia
     &                             + doziayic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + domegadxia*dzic + domegadzic*dxia
     &                             + doxiazic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + domegadyia*dzic + domegadzic*dyia
     &                             + doyiazic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + domegadzia*dzic + domegadzic*dzia
     &                             + doziazic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + domegadxia*dxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + domegadyia*dxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + domegadzia*dxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + domegadxia*dyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + domegadyia*dyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + domegadzia*dyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + domegadxia*dzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + domegadyia*dzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + domegadzia*dzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + 2.0d0 * domegadxib * dxib
     &                             + doxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + doxibyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + doxibzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + doxibyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + 2.0d0 * domegadyib * dyib
     &                             + doyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + doyibzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + doxibzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + doyibzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + 2.0 * domegadzib * dzib
     &                             + dozibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + domegadxia*dxib + domegadxib*dxia
     &                             + doxibxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + domegadxia*dyib + domegadyib*dxia
     &                             + doyibxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + domegadxia*dzib + domegadzib*dxia
     &                             + dozibxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + domegadyia*dxib + domegadxib*dyia
     &                             + doxibyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + domegadyia*dyib + domegadyib*dyia
     &                             + doyibyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + domegadyia*dzib + domegadzib*dyia
     &                             + dozibyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + domegadzia*dxib + domegadxib*dzia
     &                             + doxibzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + domegadzia*dyib + domegadyib*dzia
     &                             + doyibzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + domegadzia*dzib + domegadzib*dzia
     &                             + dozibzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + doxibxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + doyibxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + dozibxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + doxibyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + doyibyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + dozibyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + doxibzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + doyibzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + dozibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + domegadxib*dxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + domegadyib*dxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + domegadzib*dxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + domegadxib*dyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + domegadyib*dyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + domegadzib*dyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + domegadxib*dzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + domegadyib*dzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + domegadzib*dzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0 * domegadxic * dxic
     &                             + doxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + doxicyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + doxiczic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + doxicyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0 * domegadyic * dyic
     &                             + doyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + doyiczic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + doxiczic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + doyiczic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0 * domegadzic * dzic
     &                             + doziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + domegadxia*dxic + domegadxic*dxia
     &                             + doxiaxic
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + domegadxia*dyic + domegadyic*dxia
     &                             + doxiayic
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + domegadxia*dzic + domegadzic*dxia
     &                             + doxiazic
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + domegadyia*dxic + domegadxic*dyia
     &                             + doyiaxic
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + domegadyia*dyic + domegadyic*dyia
     &                             + doyiayic
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + domegadyia*dzic + domegadzic*dyia
     &                             + doyiazic
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + domegadzia*dxic + domegadxic*dzia
     &                             + doziaxic
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + domegadzia*dyic + domegadyic*dzia
     &                             + doziayic
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + domegadzia*dzic + domegadzic*dzia
     &                             + doziazic
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + doxibxic
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + doxibyic
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + doxibzic
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + doyibxic
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + doyibyic
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + doyibzic
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + dozibxic
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + dozibyic
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + dozibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + domegadxic*dxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + domegadyic*dxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + domegadzic*dxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + domegadxic*dyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + domegadyic*dyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + domegadzic*dyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + domegadxic*dzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + domegadyic*dzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + domegadzic*dzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + domegadxia*dxid
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + domegadxia*dyid
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + domegadxia*dzid
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + domegadyia*dxid
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + domegadyia*dyid
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + domegadyia*dzid
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + domegadzia*dxid
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + domegadzia*dyid
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + domegadzia*dzid
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + domegadxib*dxid
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + domegadxib*dyid
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + domegadxib*dzid
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + domegadyib*dxid
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + domegadyib*dyid
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + domegadyib*dzid
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + domegadzib*dxid
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + domegadzib*dyid
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + domegadzib*dzid
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + domegadxic*dxid
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + domegadxic*dyid
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + domegadxic*dzid
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + domegadyic*dxid
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + domegadyic*dyid
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + domegadyic*dzid
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + domegadzic*dxid
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + domegadzic*dyid
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + domegadzic*dzid
               end if
c
c     get the angle-torsion values for the second angle
c
               v1 = kant(4,iangtor)
               v2 = kant(5,iangtor)
               v3 = kant(6,iangtor)
               k = iat(3,iangtor)
               dot = xbc*xdc + ybc*ydc + zbc*zdc
               cosang = dot / sqrt(rcb2*rdc2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               dedphi = atorunit * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = atorunit * dt
     &                       * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               d2dt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dt = d2dt * fgrp
               end if
c
c     first and second derivative components for the second angle
c
               termb = -1.0d0 / (rcb2*sqrt(ru2))
               termd = 1.0d0 / (rdc2*sqrt(ru2))
               domegadxib = termb * (ybc*zu - zbc*yu)
               domegadyib = termb * (zbc*xu - xbc*zu)
               domegadzib = termb * (xbc*yu - ybc*xu)
               domegadxid = termd * (ydc*zu - zdc*yu)
               domegadyid = termd * (zdc*xu - xdc*zu)
               domegadzid = termd * (xdc*yu - ydc*xu)
               domegadxic = -domegadxib - domegadxid
               domegadyic = -domegadyib - domegadyid
               domegadzic = -domegadzib - domegadzid
c
c     abbreviations used in defining chain rule terms
c
               xrbc = 2.0d0 * xbc / rcb2
               yrbc = 2.0d0 * ybc / rcb2
               zrbc = 2.0d0 * zbc / rcb2
               xrdc = 2.0d0 * xdc / rdc2
               yrdc = 2.0d0 * ydc / rdc2
               zrdc = 2.0d0 * zdc / rdc2
               xbcp = (ybc*zu-zbc*yu) / ru2
               ybcp = (zbc*xu-xbc*zu) / ru2
               zbcp = (xbc*yu-ybc*xu) / ru2
               xdcp = (ydc*zu-zdc*yu) / ru2
               ydcp = (zdc*xu-xdc*zu) / ru2
               zdcp = (xdc*yu-ydc*xu) / ru2
c
c     chain rule terms for second derivative components
c
               doxibxib = termb*(xbc*xdc-dot) + domegadxib*(xdcp-xrbc)
               doxibyib = termb*(zu+ybc*xdc) + domegadxib*(ydcp-yrbc)
               doxibzib = termb*(zbc*xdc-yu) + domegadxib*(zdcp-zrbc)
               doyibyib = termb*(ybc*ydc-dot) + domegadyib*(ydcp-yrbc)
               doyibzib = termb*(xu+zbc*ydc) + domegadyib*(zdcp-zrbc)
               dozibzib = termb*(zbc*zdc-dot) + domegadzib*(zdcp-zrbc)
               doxidxid = termd*(dot-xbc*xdc) - domegadxid*(xbcp+xrdc)
               doxidyid = termd*(zu-ydc*xbc) - domegadxid*(ybcp+yrdc)
               doxidzid = -termd*(yu+zdc*xbc) - domegadxid*(zbcp+zrdc)
               doyidyid = termd*(dot-ybc*ydc) - domegadyid*(ybcp+yrdc)
               doyidzid = termd*(xu-zdc*ybc) - domegadyid*(zbcp+zrdc)
               dozidzid = termd*(dot-zbc*zdc) - domegadzid*(zbcp+zrdc)
               doxibxid = termb*(ybc*ybc+zbc*zbc) - domegadxib*xbcp
               doxibyid = -termb*xbc*ybc - domegadxib*ybcp
               doxibzid = -termb*xbc*zbc - domegadxib*zbcp
               doyibxid = -termb*xbc*ybc - domegadyib*xbcp
               doyibyid = termb*(xbc*xbc+zbc*zbc) - domegadyib*ybcp
               doyibzid = -termb*ybc*zbc - domegadyib*zbcp
               dozibxid = -termb*xbc*zbc - domegadzib*xbcp
               dozibyid = -termb*ybc*zbc - domegadzib*ybcp
               dozibzid = termb*(xbc*xbc+ybc*ybc) - domegadzib*zbcp
c
c     get some second derivative chain rule terms by difference
c
               doxicxib = -doxibxib - doxibxid
               doxicyib = -doxibyib - doyibxid
               doxiczib = -doxibzib - dozibxid
               doyicxib = -doxibyib - doxibyid
               doyicyib = -doyibyib - doyibyid
               doyiczib = -doyibzib - dozibyid
               dozicxib = -doxibzib - doxibzid
               dozicyib = -doyibzib - doyibzid
               doziczib = -dozibzib - dozibzid
               doxicxid = -doxidxid - doxibxid
               doxicyid = -doxidyid - doxibyid
               doxiczid = -doxidzid - doxibzid
               doyicxid = -doxidyid - doyibxid
               doyicyid = -doyidyid - doyibyid
               doyiczid = -doyidzid - doyibzid
               dozicxid = -doxidzid - dozibxid
               dozicyid = -doyidzid - dozibyid
               doziczid = -dozidzid - dozibzid
               doxicxic = -doxicxib - doxicxid
               doxicyic = -doxicyib - doxicyid
               doxiczic = -doxiczib - doxiczid
               doyicyic = -doyicyib - doyicyid
               doyiczic = -doyiczib - doyiczid
               doziczic = -doziczib - doziczid
c
c     scale the first-derivatives of the second angle
c
               domegadxib = domegadxib * radian
               domegadyib = domegadyib * radian
               domegadzib = domegadzib * radian
               domegadxid = domegadxid * radian
               domegadyid = domegadyid * radian
               domegadzid = domegadzid * radian
               domegadxic = domegadxic * radian
               domegadyic = domegadyic * radian
               domegadzic = domegadzic * radian
c
c     scale the second-derivatives of the second angle
c
               doxibxib = doxibxib * d2dt
               doxibyib = doxibyib * d2dt
               doxibzib = doxibzib * d2dt
               doyibyib = doyibyib * d2dt
               doyibzib = doyibzib * d2dt
               dozibzib = dozibzib * d2dt
               doxidxid = doxidxid * d2dt
               doxidyid = doxidyid * d2dt
               doxidzid = doxidzid * d2dt
               doyidyid = doyidyid * d2dt
               doyidzid = doyidzid * d2dt
               dozidzid = dozidzid * d2dt
               doxibxid = doxibxid * d2dt
               doxibyid = doxibyid * d2dt
               doxibzid = doxibzid * d2dt
               doyibxid = doyibxid * d2dt
               doyibyid = doyibyid * d2dt
               doyibzid = doyibzid * d2dt
               dozibxid = dozibxid * d2dt
               dozibyid = dozibyid * d2dt
               dozibzid = dozibzid * d2dt
               doxicxib = doxicxib * d2dt
               doxicyib = doxicyib * d2dt
               doxiczib = doxiczib * d2dt
               doyicxib = doyicxib * d2dt
               doyicyib = doyicyib * d2dt
               doyiczib = doyiczib * d2dt
               dozicxib = dozicxib * d2dt
               dozicyib = dozicyib * d2dt
               doziczib = doziczib * d2dt
               doxicxid = doxicxid * d2dt
               doxicyid = doxicyid * d2dt
               doxiczid = doxiczid * d2dt
               doyicxid = doyicxid * d2dt
               doyicyid = doyicyid * d2dt
               doyiczid = doyiczid * d2dt
               dozicxid = dozicxid * d2dt
               dozicyid = dozicyid * d2dt
               doziczid = doziczid * d2dt
               doxicxic = doxicxic * d2dt
               doxicyic = doxicyic * d2dt
               doxiczic = doxiczic * d2dt
               doyicyic = doyicyic * d2dt
               doyiczic = doyiczic * d2dt
               doziczic = doziczic * d2dt
c
c     chain rule terms for first derivative components
c
               dxia = dedphi * dphidxia
               dyia = dedphi * dphidyia
               dzia = dedphi * dphidzia
               dxib = dedphi * dphidxib
               dyib = dedphi * dphidyib
               dzib = dedphi * dphidzib
               dxic = dedphi * dphidxic
               dyic = dedphi * dphidyic
               dzic = dedphi * dphidzic
               dxid = dedphi * dphidxid
               dyid = dedphi * dphidyid
               dzid = dedphi * dphidzid
               dedphi = dedphi * dt
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + domegadxib*dxia
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + domegadxib*dyia
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + domegadxib*dzia
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + domegadyib*dxia
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + domegadyib*dyia
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + domegadyib*dzia
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + domegadzib*dxia
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + domegadzib*dyia
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + domegadzib*dzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + domegadxic*dxia
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + domegadxic*dyia
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + domegadxic*dzia
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + domegadyic*dxia
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + domegadyic*dyia
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + domegadyic*dzia
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + domegadzic*dxia
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + domegadzic*dyia
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + domegadzic*dzia
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + domegadxid*dxia
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + domegadxid*dyia
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + domegadxid*dzia
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + domegadyid*dxia
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + domegadyid*dyia
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + domegadyid*dzia
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + domegadzid*dxia
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + domegadzid*dyia
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + domegadzid*dzia
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + 2.0d0 * domegadxib * dxib
     &                             + doxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + doxibyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + doxibzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + doxibyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + 2.0d0 * domegadyib * dyib
     &                             + doyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + doyibzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + doxibzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + doyibzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + 2.0d0 * domegadzib * dzib
     &                             + dozibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + domegadxib*dxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + domegadyib*dxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + domegadzib*dxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + domegadxib*dyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + domegadyib*dyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + domegadzib*dyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + domegadxib*dzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + domegadyib*dzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + domegadzib*dzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + doxicxib
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + doxicyib
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + doxiczib
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + doyicxib
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + doyicyib
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + doyiczib
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + dozicxib
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + dozicyib
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + doziczib
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + domegadxib*dxid + domegadxid*dxib
     &                             + doxibxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + domegadyib*dxid + domegadxid*dyib
     &                             + doyibxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + domegadzib*dxid + domegadxid*dzib
     &                             + dozibxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + domegadxib*dyid + domegadyid*dxib
     &                             + doxibyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + domegadyib*dyid + domegadyid*dyib
     &                             + doyibyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + domegadzib*dyid + domegadyid*dzib
     &                             + dozibyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + domegadxib*dzid + domegadzid*dxib
     &                             + doxibzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + domegadyib*dzid + domegadzid*dyib
     &                             + doyibzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + domegadzib*dzid + domegadzid*dzib
     &                             + dozibzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0 * domegadxic * dxic
     &                             + doxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + doxicyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + doxiczic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + doxicyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0 * domegadyic * dyic
     &                             + doyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + doyiczic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + doxiczic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + doyiczic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0 * domegadzic * dzic
     &                             + doziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + domegadxic*dxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + domegadyic*dxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + domegadzic*dxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + domegadxic*dyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + domegadyic*dyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + domegadzic*dyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + domegadxic*dzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + domegadyic*dzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + domegadzic*dzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + doxicxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + doyicxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + dozicxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + doxicyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + doyicyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + dozicyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + doxiczib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + doyiczib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + doziczib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + domegadxic*dxid + domegadxid*dxic
     &                             + doxicxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + domegadyic*dxid + domegadxid*dyic
     &                             + doyicxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + domegadzic*dxid + domegadxid*dzic
     &                             + dozicxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + domegadxic*dyid + domegadyid*dxic
     &                             + doxicyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + domegadyic*dyid + domegadyid*dyic
     &                             + doyicyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + domegadzic*dyid + domegadyid*dzic
     &                             + dozicyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + domegadxic*dzid + domegadzid*dxic
     &                             + doxiczid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + domegadyic*dzid + domegadzid*dyic
     &                             + doyiczid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + domegadzic*dzid + domegadzid*dzic
     &                             + doziczid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + 2.0d0 * domegadxid * dxid
     &                             + doxidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + domegadxid*dyid + domegadyid*dxid
     &                             + doxidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + domegadxid*dzid + domegadzid*dxid
     &                             + doxidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + domegadxid*dyid + domegadyid*dxid
     &                             + doxidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + 2.0d0 * domegadyid * dyid
     &                             + doyidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + domegadyid*dzid + domegadzid*dyid
     &                             + doyidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + domegadxid*dzid + domegadzid*dxid
     &                             + doxidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + domegadyid*dzid + domegadzid*dyid
     &                             + doyidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + 2.0d0 * domegadzid * dzid
     &                             + dozidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + domegadxid*dxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + domegadyid*dxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + domegadzid*dxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + domegadxid*dyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + domegadyid*dyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + domegadzid*dyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + domegadxid*dzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + domegadyid*dzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + domegadzid*dzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + domegadxib*dxid + domegadxid*dxib
     &                             + doxibxid
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + domegadxib*dyid + domegadyid*dxib
     &                             + doxibyid
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + domegadxib*dzid + domegadzid*dxib
     &                             + doxibzid
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + domegadyib*dxid + domegadxid*dyib
     &                             + doyibxid
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + domegadyib*dyid + domegadyid*dyib
     &                             + doyibyid
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + domegadyib*dzid + domegadzid*dyib
     &                             + doyibzid
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + domegadzib*dxid + domegadxid*dzib
     &                             + dozibxid
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + domegadzib*dyid + domegadyid*dzib
     &                             + dozibyid
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + domegadzib*dzid + domegadzid*dzib
     &                             + dozibzid
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + domegadxic*dxid + domegadxid*dxic
     &                             + doxicxid
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + domegadxic*dyid + domegadyid*dxic
     &                             + doxicyid
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + domegadxic*dzid + domegadzid*dxic
     &                             + doxiczid
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + domegadyic*dxid + domegadxid*dyic
     &                             + doyicxid
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + domegadyic*dyid + domegadyid*dyic
     &                             + doyicyid
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + domegadyic*dzid + domegadzid*dyic
     &                             + doyiczid
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + domegadzic*dxid + domegadxid*dzic
     &                             + dozicxid
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + domegadzic*dyid + domegadyid*dzic
     &                             + dozicyid
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + domegadzic*dzid + domegadzid*dzic
     &                             + doziczid
               end if
            end if
         end if
      end do
      return
      end
