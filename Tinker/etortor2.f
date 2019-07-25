c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine etortor2  --  atomwise torsion-torsion Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "etortor2" calculates the torsion-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine etortor2 (i)
      use sizes
      use atoms
      use bitor
      use bound
      use group
      use hessn
      use ktrtor
      use math
      use torpot
      use tortor
      use units
      implicit none
      integer i,j,k,itortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer nlo,nhi,nt,xlo,ylo
      real*8 e,fgrp,sign
      real*8 angle1,angle2
      real*8 value1,value2
      real*8 cosine1,cosine2
      real*8 dedang1,dedang2
      real*8 d2eda1a1,d2eda2a2
      real*8 d2eda1a2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xed,yed,zed
      real*8 xec,yec,zec
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xv,yv,zv
      real*8 xtu,ytu,ztu
      real*8 xuv,yuv,zuv
      real*8 xh,yh,x1l,x1u,y1l,y1u
      real*8 rt2,ru2,rv2
      real*8 rtru,rcb,rurv,rdc
      real*8 da1dxt,da1dyt,da1dzt
      real*8 da1dxu,da1dyu,da1dzu
      real*8 da1dxia,da1dyia,da1dzia
      real*8 da1dxib,da1dyib,da1dzib
      real*8 da1dxic,da1dyic,da1dzic
      real*8 da1dxid,da1dyid,da1dzid
      real*8 da2dxv,da2dyv,da2dzv
      real*8 da2dxu,da2dyu,da2dzu
      real*8 da2dxie,da2dyie,da2dzie
      real*8 da2dxib,da2dyib,da2dzib
      real*8 da2dxic,da2dyic,da2dzic
      real*8 da2dxid,da2dyid,da2dzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 xydc2,xzdc2,yzdc2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 rdcxv,rdcyv,rdczv,rdcv2
      real*8 rdcxu,rdcyu,rdczu,rdcu2
      real*8 da1dxibt,da1dyibt,da1dzibt
      real*8 da1dxibu,da1dyibu,da1dzibu
      real*8 da1dxict,da1dyict,da1dzict
      real*8 da1dxicu,da1dyicu,da1dzicu
      real*8 da2dxidu,da2dyidu,da2dzidu
      real*8 da2dxidv,da2dyidv,da2dzidv
      real*8 da2dxicu,da2dyicu,da2dzicu
      real*8 da2dxicv,da2dyicv,da2dzicv
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
      real*8 dxibxib2,dyibyib2,dzibzib2
      real*8 dxicxic2,dyicyic2,dziczic2
      real*8 dxidxid2,dyidyid2,dzidzid2
      real*8 dxiexie2,dyieyie2,dziezie2
      real*8 dxibyib2,dxibzib2,dyibzib2
      real*8 dxicyic2,dxiczic2,dyiczic2
      real*8 dxidyid2,dxidzid2,dyidzid2
      real*8 dxieyie2,dxiezie2,dyiezie2
      real*8 dxibxic2,dxibyic2,dxibzic2
      real*8 dyibxic2,dyibyic2,dyibzic2
      real*8 dzibxic2,dzibyic2,dzibzic2
      real*8 dxibxid2,dxibyid2,dxibzid2
      real*8 dyibxid2,dyibyid2,dyibzid2
      real*8 dzibxid2,dzibyid2,dzibzid2
      real*8 dxibxie2,dxibyie2,dxibzie2
      real*8 dyibxie2,dyibyie2,dyibzie2
      real*8 dzibxie2,dzibyie2,dzibzie2
      real*8 dxicxid2,dxicyid2,dxiczid2
      real*8 dyicxid2,dyicyid2,dyiczid2
      real*8 dzicxid2,dzicyid2,dziczid2
      real*8 dxicxie2,dxicyie2,dxiczie2
      real*8 dyicxie2,dyicyie2,dyiczie2
      real*8 dzicxie2,dzicyie2,dziczie2
      real*8 dxidxie2,dxidyie2,dxidzie2
      real*8 dyidxie2,dyidyie2,dyidzie2
      real*8 dzidxie2,dzidyie2,dzidzie2
      real*8 ftt(4),ft12(4)
      real*8 ft1(4),ft2(4)
      logical proceed
c
c
c     compute the Hessian elements of the torsion-torsions
c
      do itortor = 1, ntortor
         j = itt(1,itortor)
         k = itt(2,itortor)
         if (itt(3,itortor) .eq. 1) then
            ia = ibitor(1,j)
            ib = ibitor(2,j)
            ic = ibitor(3,j)
            id = ibitor(4,j)
            ie = ibitor(5,j)
         else
            ia = ibitor(5,j)
            ib = ibitor(4,j)
            ic = ibitor(3,j)
            id = ibitor(2,j)
            ie = ibitor(1,j)
         end if
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic
     &                 .or. i.eq.id .or. i.eq.ie)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
c
c     compute the values of the torsional angles
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
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xed,yed,zed)
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
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            xuv = yu*zv - yv*zu
            yuv = zu*xv - zv*xu
            zuv = xu*yv - xv*yu
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru .ne. 0.0d0 .and. rurv.ne.0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0d0,max(-1.0d0,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0d0)  angle1 = -angle1
               value1 = angle1
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0d0,max(-1.0d0,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0d0)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               ylo = nlo
               xh = ttx(xlo+1,k) - ttx(xlo,k)
               yh = tty(ylo+1,k) - tty(ylo,k)
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint2 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,
     &                       value1,value2,e,dedang1,dedang2,
     &                       d2eda1a2,d2eda1a1,d2eda2a2)
               dedang1 = sign * ttorunit * radian * dedang1
               dedang2 = sign * ttorunit * radian * dedang2
               d2eda1a1 = ttorunit * radian**2 * d2eda1a1
               d2eda2a2 = ttorunit * radian**2 * d2eda2a2
               d2eda1a2 = ttorunit * radian**2 * d2eda1a2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedang1 = dedang1 * fgrp
                  dedang2 = dedang2 * fgrp
                  d2eda1a1 = d2eda1a1 * fgrp
                  d2eda2a2 = d2eda2a2 * fgrp
                  d2eda1a2 = d2eda1a2 * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
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
               da1dxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               da1dyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               da1dzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               da1dxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               da1dyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               da1dzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * da1dxt
               rcbyt = -2.0d0 * rcb * da1dyt
               rcbzt = -2.0d0 * rcb * da1dzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * da1dxu
               rcbyu = 2.0d0 * rcb * da1dyu
               rcbzu = 2.0d0 * rcb * da1dzu
               rcbu2 = rcb * ru2
               da1dxibt = yca*da1dzt - zca*da1dyt
               da1dxibu = zdc*da1dyu - ydc*da1dzu
               da1dyibt = zca*da1dxt - xca*da1dzt
               da1dyibu = xdc*da1dzu - zdc*da1dxu
               da1dzibt = xca*da1dyt - yca*da1dxt
               da1dzibu = ydc*da1dxu - xdc*da1dyu
               da1dxict = zba*da1dyt - yba*da1dzt
               da1dxicu = ydb*da1dzu - zdb*da1dyu
               da1dyict = xba*da1dzt - zba*da1dxt
               da1dyicu = zdb*da1dxu - xdb*da1dzu
               da1dzict = yba*da1dxt - xba*da1dyt
               da1dzicu = xdb*da1dyu - ydb*da1dxu
c
c     chain rule terms for first derivative components
c
               da1dxia = zcb*da1dyt - ycb*da1dzt
               da1dyia = xcb*da1dzt - zcb*da1dxt
               da1dzia = ycb*da1dxt - xcb*da1dyt
               da1dxib = da1dxibt + da1dxibu
               da1dyib = da1dyibt + da1dyibu
               da1dzib = da1dzibt + da1dzibu
               da1dxic = da1dxict + da1dxicu
               da1dyic = da1dyict + da1dyicu
               da1dzic = da1dzict + da1dzicu
               da1dxid = zcb*da1dyu - ycb*da1dzu
               da1dyid = xcb*da1dzu - zcb*da1dxu
               da1dzid = ycb*da1dxu - xcb*da1dyu
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*da1dxia
               dxiayia = rcbxt*da1dyia - zcb*rcb/rt2
               dxiazia = rcbxt*da1dzia + ycb*rcb/rt2
               dxiaxic = rcbxt*da1dxict + xcb*xt/rcbt2
               dxiayic = rcbxt*da1dyict - da1dzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*da1dzict + da1dyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*da1dyia
               dyiazia = rcbyt*da1dzia - xcb*rcb/rt2
               dyiaxib = rcbyt*da1dxibt - da1dzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiaxic = rcbyt*da1dxict + da1dzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*da1dyict + ycb*yt/rcbt2
               dyiazic = rcbyt*da1dzict - da1dxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*da1dzia
               dziaxib = rcbzt*da1dxibt + da1dyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*da1dyibt - da1dxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziaxic = rcbzt*da1dxict - da1dyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*da1dyict + da1dxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*da1dzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*da1dxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*da1dxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*da1dxibu/ru2
               dxibyic = -ycb*da1dxib/(rcb*rcb) + da1dzt + da1dzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*da1dxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*da1dxibu/ru2
               dxibxid = rcbxu*da1dxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*da1dxibu - da1dzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*da1dxibu + da1dyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*da1dzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*da1dzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*da1dzibu/ru2
               dyibxic = -xcb*da1dyib/(rcb*rcb) - da1dzt - da1dzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*da1dyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*da1dyibu/ru2
               dyibyic = -ycb*da1dyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*da1dyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*da1dyibu/ru2
               dyibxid = rcbxu*da1dyibu + da1dzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*da1dyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*da1dyibu - da1dxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*da1dzib/(rcb*rcb) + da1dyt + da1dyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*da1dzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*da1dzibu/ru2
               dzibzic = -zcb*da1dzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*da1dzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*da1dzibu/ru2
               dzibxid = rcbxu*da1dzibu - da1dyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*da1dzibu + da1dxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*da1dzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*da1dxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*da1dxicu + da1dzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*da1dxicu - da1dyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*da1dyicu - da1dzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*da1dyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*da1dyicu + da1dxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*da1dzicu + da1dyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*da1dzicu - da1dxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*da1dzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*da1dxid
               dxidyid = rcbxu*da1dyid + zcb*rcb/ru2
               dxidzid = rcbxu*da1dzid - ycb*rcb/ru2
               dyidyid = rcbyu*da1dyid
               dyidzid = rcbyu*da1dzid + xcb*rcb/ru2
               dzidzid = rcbzu*da1dzid
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
                  hessx(1,ia) = hessx(1,ia) + dedang1*dxiaxia
     &                             + d2eda1a1*da1dxia*da1dxia
                  hessy(1,ia) = hessy(1,ia) + dedang1*dxiayia
     &                             + d2eda1a1*da1dxia*da1dyia
                  hessz(1,ia) = hessz(1,ia) + dedang1*dxiazia
     &                             + d2eda1a1*da1dxia*da1dzia
                  hessx(2,ia) = hessx(2,ia) + dedang1*dxiayia
     &                             + d2eda1a1*da1dxia*da1dyia
                  hessy(2,ia) = hessy(2,ia) + dedang1*dyiayia
     &                             + d2eda1a1*da1dyia*da1dyia
                  hessz(2,ia) = hessz(2,ia) + dedang1*dyiazia
     &                             + d2eda1a1*da1dyia*da1dzia
                  hessx(3,ia) = hessx(3,ia) + dedang1*dxiazia
     &                             + d2eda1a1*da1dxia*da1dzia
                  hessy(3,ia) = hessy(3,ia) + dedang1*dyiazia
     &                             + d2eda1a1*da1dyia*da1dzia
                  hessz(3,ia) = hessz(3,ia) + dedang1*dziazia
     &                             + d2eda1a1*da1dzia*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang1*dxiaxib
     &                             + d2eda1a1*da1dxia*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang1*dyiaxib
     &                             + d2eda1a1*da1dyia*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang1*dziaxib
     &                             + d2eda1a1*da1dzia*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang1*dxiayib
     &                             + d2eda1a1*da1dxia*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang1*dyiayib
     &                             + d2eda1a1*da1dyia*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang1*dziayib
     &                             + d2eda1a1*da1dzia*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang1*dxiazib
     &                             + d2eda1a1*da1dxia*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang1*dyiazib
     &                             + d2eda1a1*da1dyia*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang1*dziazib
     &                             + d2eda1a1*da1dzia*da1dzib
                  hessx(1,ic) = hessx(1,ic) + dedang1*dxiaxic
     &                             + d2eda1a1*da1dxia*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang1*dyiaxic
     &                             + d2eda1a1*da1dyia*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang1*dziaxic
     &                             + d2eda1a1*da1dzia*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang1*dxiayic
     &                             + d2eda1a1*da1dxia*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang1*dyiayic
     &                             + d2eda1a1*da1dyia*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang1*dziayic
     &                             + d2eda1a1*da1dzia*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang1*dxiazic
     &                             + d2eda1a1*da1dxia*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang1*dyiazic
     &                             + d2eda1a1*da1dyia*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang1*dziazic
     &                             + d2eda1a1*da1dzia*da1dzic
                  hessx(1,id) = hessx(1,id) + dedang1*dxiaxid
     &                             + d2eda1a1*da1dxia*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang1*dyiaxid
     &                             + d2eda1a1*da1dyia*da1dxid
                  hessz(1,id) = hessz(1,id) + dedang1*dziaxid
     &                             + d2eda1a1*da1dzia*da1dxid
                  hessx(2,id) = hessx(2,id) + dedang1*dxiayid
     &                             + d2eda1a1*da1dxia*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang1*dyiayid
     &                             + d2eda1a1*da1dyia*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang1*dziayid
     &                             + d2eda1a1*da1dzia*da1dyid
                  hessx(3,id) = hessx(3,id) + dedang1*dxiazid
     &                             + d2eda1a1*da1dxia*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang1*dyiazid
     &                             + d2eda1a1*da1dyia*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang1*dziazid
     &                             + d2eda1a1*da1dzia*da1dzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedang1*dxibxib
     &                             + d2eda1a1*da1dxib*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang1*dxibyib
     &                             + d2eda1a1*da1dxib*da1dyib
                  hessz(1,ib) = hessz(1,ib) + dedang1*dxibzib
     &                             + d2eda1a1*da1dxib*da1dzib
                  hessx(2,ib) = hessx(2,ib) + dedang1*dxibyib
     &                             + d2eda1a1*da1dxib*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang1*dyibyib
     &                             + d2eda1a1*da1dyib*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang1*dyibzib
     &                             + d2eda1a1*da1dyib*da1dzib
                  hessx(3,ib) = hessx(3,ib) + dedang1*dxibzib
     &                             + d2eda1a1*da1dxib*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang1*dyibzib
     &                             + d2eda1a1*da1dyib*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang1*dzibzib
     &                             + d2eda1a1*da1dzib*da1dzib
                  hessx(1,ia) = hessx(1,ia) + dedang1*dxiaxib
     &                             + d2eda1a1*da1dxib*da1dxia
                  hessy(1,ia) = hessy(1,ia) + dedang1*dxiayib
     &                             + d2eda1a1*da1dyib*da1dxia
                  hessz(1,ia) = hessz(1,ia) + dedang1*dxiazib
     &                             + d2eda1a1*da1dzib*da1dxia
                  hessx(2,ia) = hessx(2,ia) + dedang1*dyiaxib
     &                             + d2eda1a1*da1dxib*da1dyia
                  hessy(2,ia) = hessy(2,ia) + dedang1*dyiayib
     &                             + d2eda1a1*da1dyib*da1dyia
                  hessz(2,ia) = hessz(2,ia) + dedang1*dyiazib
     &                             + d2eda1a1*da1dzib*da1dyia
                  hessx(3,ia) = hessx(3,ia) + dedang1*dziaxib
     &                             + d2eda1a1*da1dxib*da1dzia
                  hessy(3,ia) = hessy(3,ia) + dedang1*dziayib
     &                             + d2eda1a1*da1dyib*da1dzia
                  hessz(3,ia) = hessz(3,ia) + dedang1*dziazib
     &                             + d2eda1a1*da1dzib*da1dzia
                  hessx(1,ic) = hessx(1,ic) + dedang1*dxibxic
     &                             + d2eda1a1*da1dxib*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang1*dyibxic
     &                             + d2eda1a1*da1dyib*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang1*dzibxic
     &                             + d2eda1a1*da1dzib*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang1*dxibyic
     &                             + d2eda1a1*da1dxib*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang1*dyibyic
     &                             + d2eda1a1*da1dyib*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang1*dzibyic
     &                             + d2eda1a1*da1dzib*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang1*dxibzic
     &                             + d2eda1a1*da1dxib*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang1*dyibzic
     &                             + d2eda1a1*da1dyib*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang1*dzibzic
     &                             + d2eda1a1*da1dzib*da1dzic
                  hessx(1,id) = hessx(1,id) + dedang1*dxibxid
     &                             + d2eda1a1*da1dxib*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang1*dyibxid
     &                             + d2eda1a1*da1dyib*da1dxid
                  hessz(1,id) = hessz(1,id) + dedang1*dzibxid
     &                             + d2eda1a1*da1dzib*da1dxid
                  hessx(2,id) = hessx(2,id) + dedang1*dxibyid
     &                             + d2eda1a1*da1dxib*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang1*dyibyid
     &                             + d2eda1a1*da1dyib*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang1*dzibyid
     &                             + d2eda1a1*da1dzib*da1dyid
                  hessx(3,id) = hessx(3,id) + dedang1*dxibzid
     &                             + d2eda1a1*da1dxib*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang1*dyibzid
     &                             + d2eda1a1*da1dyib*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang1*dzibzid
     &                             + d2eda1a1*da1dzib*da1dzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedang1*dxicxic
     &                             + d2eda1a1*da1dxic*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang1*dxicyic
     &                             + d2eda1a1*da1dxic*da1dyic
                  hessz(1,ic) = hessz(1,ic) + dedang1*dxiczic
     &                             + d2eda1a1*da1dxic*da1dzic
                  hessx(2,ic) = hessx(2,ic) + dedang1*dxicyic
     &                             + d2eda1a1*da1dxic*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang1*dyicyic
     &                             + d2eda1a1*da1dyic*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang1*dyiczic
     &                             + d2eda1a1*da1dyic*da1dzic
                  hessx(3,ic) = hessx(3,ic) + dedang1*dxiczic
     &                             + d2eda1a1*da1dxic*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang1*dyiczic
     &                             + d2eda1a1*da1dyic*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang1*dziczic
     &                             + d2eda1a1*da1dzic*da1dzic
                  hessx(1,ia) = hessx(1,ia) + dedang1*dxiaxic
     &                             + d2eda1a1*da1dxic*da1dxia
                  hessy(1,ia) = hessy(1,ia) + dedang1*dxiayic
     &                             + d2eda1a1*da1dyic*da1dxia
                  hessz(1,ia) = hessz(1,ia) + dedang1*dxiazic
     &                             + d2eda1a1*da1dzic*da1dxia
                  hessx(2,ia) = hessx(2,ia) + dedang1*dyiaxic
     &                             + d2eda1a1*da1dxic*da1dyia
                  hessy(2,ia) = hessy(2,ia) + dedang1*dyiayic
     &                             + d2eda1a1*da1dyic*da1dyia
                  hessz(2,ia) = hessz(2,ia) + dedang1*dyiazic
     &                             + d2eda1a1*da1dzic*da1dyia
                  hessx(3,ia) = hessx(3,ia) + dedang1*dziaxic
     &                             + d2eda1a1*da1dxic*da1dzia
                  hessy(3,ia) = hessy(3,ia) + dedang1*dziayic
     &                             + d2eda1a1*da1dyic*da1dzia
                  hessz(3,ia) = hessz(3,ia) + dedang1*dziazic
     &                             + d2eda1a1*da1dzic*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang1*dxibxic
     &                             + d2eda1a1*da1dxic*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang1*dxibyic
     &                             + d2eda1a1*da1dyic*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang1*dxibzic
     &                             + d2eda1a1*da1dzic*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang1*dyibxic
     &                             + d2eda1a1*da1dxic*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang1*dyibyic
     &                             + d2eda1a1*da1dyic*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang1*dyibzic
     &                             + d2eda1a1*da1dzic*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang1*dzibxic
     &                             + d2eda1a1*da1dxic*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang1*dzibyic
     &                             + d2eda1a1*da1dyic*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang1*dzibzic
     &                             + d2eda1a1*da1dzic*da1dzib
                  hessx(1,id) = hessx(1,id) + dedang1*dxicxid
     &                             + d2eda1a1*da1dxic*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang1*dyicxid
     &                             + d2eda1a1*da1dyic*da1dxid
                  hessz(1,id) = hessz(1,id) + dedang1*dzicxid
     &                             + d2eda1a1*da1dzic*da1dxid
                  hessx(2,id) = hessx(2,id) + dedang1*dxicyid
     &                             + d2eda1a1*da1dxic*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang1*dyicyid
     &                             + d2eda1a1*da1dyic*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang1*dzicyid
     &                             + d2eda1a1*da1dzic*da1dyid
                  hessx(3,id) = hessx(3,id) + dedang1*dxiczid
     &                             + d2eda1a1*da1dxic*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang1*dyiczid
     &                             + d2eda1a1*da1dyic*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang1*dziczid
     &                             + d2eda1a1*da1dzic*da1dzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedang1*dxidxid
     &                             + d2eda1a1*da1dxid*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang1*dxidyid
     &                             + d2eda1a1*da1dxid*da1dyid
                  hessz(1,id) = hessz(1,id) + dedang1*dxidzid
     &                             + d2eda1a1*da1dxid*da1dzid
                  hessx(2,id) = hessx(2,id) + dedang1*dxidyid
     &                             + d2eda1a1*da1dxid*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang1*dyidyid
     &                             + d2eda1a1*da1dyid*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang1*dyidzid
     &                             + d2eda1a1*da1dyid*da1dzid
                  hessx(3,id) = hessx(3,id) + dedang1*dxidzid
     &                             + d2eda1a1*da1dxid*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang1*dyidzid
     &                             + d2eda1a1*da1dyid*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang1*dzidzid
     &                             + d2eda1a1*da1dzid*da1dzid
                  hessx(1,ia) = hessx(1,ia) + dedang1*dxiaxid
     &                             + d2eda1a1*da1dxid*da1dxia
                  hessy(1,ia) = hessy(1,ia) + dedang1*dxiayid
     &                             + d2eda1a1*da1dyid*da1dxia
                  hessz(1,ia) = hessz(1,ia) + dedang1*dxiazid
     &                             + d2eda1a1*da1dzid*da1dxia
                  hessx(2,ia) = hessx(2,ia) + dedang1*dyiaxid
     &                             + d2eda1a1*da1dxid*da1dyia
                  hessy(2,ia) = hessy(2,ia) + dedang1*dyiayid
     &                             + d2eda1a1*da1dyid*da1dyia
                  hessz(2,ia) = hessz(2,ia) + dedang1*dyiazid
     &                             + d2eda1a1*da1dzid*da1dyia
                  hessx(3,ia) = hessx(3,ia) + dedang1*dziaxid
     &                             + d2eda1a1*da1dxid*da1dzia
                  hessy(3,ia) = hessy(3,ia) + dedang1*dziayid
     &                             + d2eda1a1*da1dyid*da1dzia
                  hessz(3,ia) = hessz(3,ia) + dedang1*dziazid
     &                             + d2eda1a1*da1dzid*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang1*dxibxid
     &                             + d2eda1a1*da1dxid*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang1*dxibyid
     &                             + d2eda1a1*da1dyid*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang1*dxibzid
     &                             + d2eda1a1*da1dzid*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang1*dyibxid
     &                             + d2eda1a1*da1dxid*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang1*dyibyid
     &                             + d2eda1a1*da1dyid*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang1*dyibzid
     &                             + d2eda1a1*da1dzid*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang1*dzibxid
     &                             + d2eda1a1*da1dxid*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang1*dzibyid
     &                             + d2eda1a1*da1dyid*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang1*dzibzid
     &                             + d2eda1a1*da1dzid*da1dzib
                  hessx(1,ic) = hessx(1,ic) + dedang1*dxicxid
     &                             + d2eda1a1*da1dxid*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang1*dxicyid
     &                             + d2eda1a1*da1dyid*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang1*dxiczid
     &                             + d2eda1a1*da1dzid*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang1*dyicxid
     &                             + d2eda1a1*da1dxid*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang1*dyicyid
     &                             + d2eda1a1*da1dyid*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang1*dyiczid
     &                             + d2eda1a1*da1dzid*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang1*dzicxid
     &                             + d2eda1a1*da1dxid*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang1*dzicyid
     &                             + d2eda1a1*da1dyid*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang1*dziczid
     &                             + d2eda1a1*da1dzid*da1dzic
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xec = xie - xic
               yec = yie - yic
               zec = zie - zic
               if (use_polymer) then
                  call image (xec,yec,zec)
               end if
               da2dxu = (yu*zdc - ydc*zu) / (ru2*rdc)
               da2dyu = (zu*xdc - zdc*xu) / (ru2*rdc)
               da2dzu = (xu*ydc - xdc*yu) / (ru2*rdc)
               da2dxv = -(yv*zdc - ydc*zv) / (rv2*rdc)
               da2dyv = -(zv*xdc - zdc*xv) / (rv2*rdc)
               da2dzv = -(xv*ydc - xdc*yv) / (rv2*rdc)
c
c     abbreviations for second derivative chain rule terms
c
               xydc2 = xdc*xdc + ydc*ydc
               xzdc2 = xdc*xdc + zdc*zdc
               yzdc2 = ydc*ydc + zdc*zdc
               rdcxu = -2.0d0 * rdc * da2dxu
               rdcyu = -2.0d0 * rdc * da2dyu
               rdczu = -2.0d0 * rdc * da2dzu
               rdcu2 = rdc * ru2
               rdcxv = 2.0d0 * rdc * da2dxv
               rdcyv = 2.0d0 * rdc * da2dyv
               rdczv = 2.0d0 * rdc * da2dzv
               rdcv2 = rdc * rv2
               da2dxicu = ydb*da2dzu - zdb*da2dyu
               da2dxicv = zed*da2dyv - yed*da2dzv
               da2dyicu = zdb*da2dxu - xdb*da2dzu
               da2dyicv = xed*da2dzv - zed*da2dxv
               da2dzicu = xdb*da2dyu - ydb*da2dxu
               da2dzicv = yed*da2dxv - xed*da2dyv
               da2dxidu = zcb*da2dyu - ycb*da2dzu
               da2dxidv = yec*da2dzv - zec*da2dyv
               da2dyidu = xcb*da2dzu - zcb*da2dxu
               da2dyidv = zec*da2dxv - xec*da2dzv
               da2dzidu = ycb*da2dxu - xcb*da2dyu
               da2dzidv = xec*da2dyv - yec*da2dxv
c
c     chain rule terms for first derivative components
c
               da2dxib = zdc*da2dyu - ydc*da2dzu
               da2dyib = xdc*da2dzu - zdc*da2dxu
               da2dzib = ydc*da2dxu - xdc*da2dyu
               da2dxic = da2dxicu + da2dxicv
               da2dyic = da2dyicu + da2dyicv
               da2dzic = da2dzicu + da2dzicv
               da2dxid = da2dxidu + da2dxidv
               da2dyid = da2dyidu + da2dyidv
               da2dzid = da2dzidu + da2dzidv
               da2dxie = zdc*da2dyv - ydc*da2dzv
               da2dyie = xdc*da2dzv - zdc*da2dxv
               da2dzie = ydc*da2dxv - xdc*da2dyv
c
c     chain rule terms for second derivative components
c
               dxibxib2 = rdcxu*da2dxib
               dxibyib2 = rdcxu*da2dyib - zdc*rdc/ru2
               dxibzib2 = rdcxu*da2dzib + ydc*rdc/ru2
               dxibxid2 = rdcxu*da2dxidu + xdc*xu/rdcu2
               dxibyid2 = rdcxu*da2dyidu - da2dzu
     &                       - (xcb*zdc*xdc+zcb*yzdc2)/rdcu2
               dxibzid2 = rdcxu*da2dzidu + da2dyu
     &                       + (xcb*ydc*xdc+ycb*yzdc2)/rdcu2
               dxibxie2 = 0.0d0
               dxibyie2 = 0.0d0
               dxibzie2 = 0.0d0
               dyibyib2 = rdcyu*da2dyib
               dyibzib2 = rdcyu*da2dzib - xdc*rdc/ru2
               dyibxic2 = rdcyu*da2dxicu - da2dzu
     &                       - (ydb*zdc*ydc+zdb*xzdc2)/rdcu2
               dyibxid2 = rdcyu*da2dxidu + da2dzu
     &                       + (ycb*zdc*ydc+zcb*xzdc2)/rdcu2
               dyibyid2 = rdcyu*da2dyidu + ydc*yu/rdcu2
               dyibzid2 = rdcyu*da2dzidu - da2dxu
     &                       - (ycb*xdc*ydc+xcb*xzdc2)/rdcu2
               dyibxie2 = 0.0d0
               dyibyie2 = 0.0d0
               dyibzie2 = 0.0d0
               dzibzib2 = rdczu*da2dzib
               dzibxic2 = rdczu*da2dxicu + da2dyu
     &                       + (zdb*ydc*zdc+ydb*xydc2)/rdcu2
               dzibyic2 = rdczu*da2dyicu - da2dxu
     &                       - (zdb*xdc*zdc+xdb*xydc2)/rdcu2
               dzibxid2 = rdczu*da2dxidu - da2dyu
     &                       - (zcb*ydc*zdc+ycb*xydc2)/rdcu2
               dzibyid2 = rdczu*da2dyidu + da2dxu
     &                       + (zcb*xdc*zdc+xcb*xydc2)/rdcu2
               dzibzid2 = rdczu*da2dzidu + zdc*zu/rdcu2
               dzibxie2 = 0.0d0
               dzibyie2 = 0.0d0
               dzibzie2 = 0.0d0
               dxicxid2 = -xdc*da2dxic/(rdc*rdc)
     &             - (ydb*(zcb*xdc+yu)-zdb*(ycb*xdc-zu))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*da2dxicu/ru2
     &             - (zed*(yec*xdc+zv)-yed*(zec*xdc-yv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*da2dxicv/rv2
               dxicyid2 = -ydc*da2dxic/(rdc*rdc) + da2dzu + da2dzv
     &             - (ydb*(zcb*ydc-xu)+zdb*(xcb*xdc+zdc*zcb))/rdcu2
     &             - 2.0d0*(zu*xcb-zcb*xu)*da2dxicu/ru2
     &             + (zed*(xec*xdc+zdc*zec)+yed*(zec*ydc+xv))/rdcv2
     &             + 2.0d0*(zv*xec-zec*xv)*da2dxicv/rv2
               dxicxie2 = rdcxv*da2dxicv + xdc*xv/rdcv2
               dxicyie2 = rdcyv*da2dxicv - da2dzv
     &                       - (yed*zdc*ydc+zed*xzdc2)/rdcv2
               dxiczie2 = rdczv*da2dxicv + da2dyv
     &                       + (zed*ydc*zdc+yed*xydc2)/rdcv2
               dyiczic2 = ydc*da2dzic/(rdc*rdc)
     &             - (xdb*(xdb*xdc+zdc*zdb)+ydb*(ydc*xdb+zu))/rdcu2
     &             - 2.0d0*(xu*zdb-xdb*zu)*da2dzicu/ru2
     &             + (yed*(xed*ydc-zv)+xed*(xed*xdc+zdc*zed))/rdcv2
     &             + 2.0d0*(xv*zed-xed*zv)*da2dzicv/rv2
               dyicxid2 = -xdc*da2dyic/(rdc*rdc) - da2dzu - da2dzv
     &             + (xdb*(zcb*xdc+yu)+zdb*(zcb*zdc+ydc*ycb))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*da2dyicu/ru2
     &             - (zed*(zec*zdc+ydc*yec)+xed*(zec*xdc-yv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*da2dyicv/rv2
               dyicyid2 = -ydc*da2dyic/(rdc*rdc)
     &             - (zdb*(xcb*ydc+zu)-xdb*(zcb*ydc-xu))/rdcu2
     &             - 2.0d0*(zu*xcb-zcb*xu)*da2dyicu/ru2
     &             - (xed*(zec*ydc+xv)-zed*(xec*ydc-zv))/rdcv2
     &             + 2.0d0*(zv*xec-zec*xv)*da2dyicv/rv2
               dyicxie2 = rdcxv*da2dyicv + da2dzv
     &                       + (xed*zdc*xdc+zed*yzdc2)/rdcv2
               dyicyie2 = rdcyv*da2dyicv + ydc*yv/rdcv2
               dyiczie2 = rdczv*da2dyicv - da2dxv
     &                       - (zed*xdc*zdc+xed*xydc2)/rdcv2
               dzicxid2 = -xdc*da2dzic/(rdc*rdc) + da2dyu + da2dyv
     &             - (xdb*(ycb*xdc-zu)+ydb*(zcb*zdc+ydc*ycb))/rdcu2
     &             - 2.0d0*(yu*zcb-ycb*zu)*da2dzicu/ru2
     &             + (yed*(zec*zdc+ydc*yec)+xed*(yec*xdc+zv))/rdcv2
     &             + 2.0d0*(yv*zec-yec*zv)*da2dzicv/rv2
               dziczid2 = -zdc*da2dzic/(rdc*rdc)
     &             - (xdb*(ycb*zdc+xu)-ydb*(xcb*zdc-yu))/rdcu2
     &             - 2.0d0*(xu*ycb-xcb*yu)*da2dzicu/ru2
     &             - (yed*(xec*zdc+yv)-xed*(yec*zdc-xv))/rdcv2
     &             + 2.0d0*(xv*yec-xec*yv)*da2dzicv/rv2
               dzicxie2 = rdcxv*da2dzicv - da2dyv
     &                       - (xed*ydc*xdc+yed*yzdc2)/rdcv2
               dzicyie2 = rdcyv*da2dzicv + da2dxv
     &                       + (yed*xdc*ydc+xed*xzdc2)/rdcv2
               dziczie2 = rdczv*da2dzicv + zdc*zv/rdcv2
               dxidxie2 = rdcxv*da2dxidv - xdc*(zec*ydc-yec*zdc)/rdcv2
               dxidyie2 = rdcyv*da2dxidv + da2dzv
     &                       + (yec*zdc*ydc+zec*xzdc2)/rdcv2
               dxidzie2 = rdczv*da2dxidv - da2dyv
     &                       - (zec*ydc*zdc+yec*xydc2)/rdcv2
               dyidxie2 = rdcxv*da2dyidv - da2dzv
     &                       - (xec*zdc*xdc+zec*yzdc2)/rdcv2
               dyidyie2 = rdcyv*da2dyidv - ydc*(xec*zdc-zec*xdc)/rdcv2
               dyidzie2 = rdczv*da2dyidv + da2dxv
     &                       + (zec*xdc*zdc+xec*xydc2)/rdcv2
               dzidxie2 = rdcxv*da2dzidv + da2dyv
     &                       + (xec*ydc*xdc+yec*yzdc2)/rdcv2
               dzidyie2 = rdcyv*da2dzidv - da2dxv
     &                       - (yec*xdc*ydc+xec*xzdc2)/rdcv2
               dzidzie2 = rdczv*da2dzidv - zdc*(yec*xdc-xec*ydc)/rdcv2
               dxiexie2 = rdcxv*da2dxie
               dxieyie2 = rdcxv*da2dyie + zdc*rdc/rv2
               dxiezie2 = rdcxv*da2dzie - ydc*rdc/rv2
               dyieyie2 = rdcyv*da2dyie
               dyiezie2 = rdcyv*da2dzie + xdc*rdc/rv2
               dziezie2 = rdczv*da2dzie
c
c     get some second derivative chain rule terms by difference
c
               dxibxic2 = -dxibxib2 - dxibxid2 - dxibxie2
               dxibyic2 = -dxibyib2 - dxibyid2 - dxibyie2
               dxibzic2 = -dxibzib2 - dxibzid2 - dxibzie2
               dyibyic2 = -dyibyib2 - dyibyid2 - dyibyie2
               dyibzic2 = -dyibzib2 - dyibzid2 - dyibzie2
               dzibzic2 = -dzibzib2 - dzibzid2 - dzibzie2
               dxicxic2 = -dxibxic2 - dxicxid2 - dxicxie2
               dxicyic2 = -dyibxic2 - dxicyid2 - dxicyie2
               dxiczic2 = -dxibzic2 - dzicxid2 - dzicxie2
               dxiczid2 = -dzibxic2 - dxiczic2 - dxiczie2
               dyicyic2 = -dyibyic2 - dyicyid2 - dyicyie2
               dyiczid2 = -dzibyic2 - dyiczic2 - dyiczie2
               dziczic2 = -dzibzic2 - dziczid2 - dziczie2
               dzicyid2 = -dyibzic2 - dyiczic2 - dzicyie2
               dxidxid2 = -dxibxid2 - dxicxid2 - dxidxie2
               dxidyid2 = -dyibxid2 - dyicxid2 - dxidyie2
               dxidzid2 = -dzibxid2 - dzicxid2 - dxidzie2
               dyidyid2 = -dyibyid2 - dyicyid2 - dyidyie2
               dyidzid2 = -dzibyid2 - dzicyid2 - dyidzie2
               dzidzid2 = -dzibzid2 - dziczid2 - dzidzie2
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ib) = hessx(1,ib) + d2eda1a2*da1dxia*da2dxib
                  hessy(1,ib) = hessy(1,ib) + d2eda1a2*da1dyia*da2dxib
                  hessz(1,ib) = hessz(1,ib) + d2eda1a2*da1dzia*da2dxib
                  hessx(2,ib) = hessx(2,ib) + d2eda1a2*da1dxia*da2dyib
                  hessy(2,ib) = hessy(2,ib) + d2eda1a2*da1dyia*da2dyib
                  hessz(2,ib) = hessz(2,ib) + d2eda1a2*da1dzia*da2dyib
                  hessx(3,ib) = hessx(3,ib) + d2eda1a2*da1dxia*da2dzib
                  hessy(3,ib) = hessy(3,ib) + d2eda1a2*da1dyia*da2dzib
                  hessz(3,ib) = hessz(3,ib) + d2eda1a2*da1dzia*da2dzib
                  hessx(1,ic) = hessx(1,ic) + d2eda1a2*da1dxia*da2dxic
                  hessy(1,ic) = hessy(1,ic) + d2eda1a2*da1dyia*da2dxic
                  hessz(1,ic) = hessz(1,ic) + d2eda1a2*da1dzia*da2dxic
                  hessx(2,ic) = hessx(2,ic) + d2eda1a2*da1dxia*da2dyic
                  hessy(2,ic) = hessy(2,ic) + d2eda1a2*da1dyia*da2dyic
                  hessz(2,ic) = hessz(2,ic) + d2eda1a2*da1dzia*da2dyic
                  hessx(3,ic) = hessx(3,ic) + d2eda1a2*da1dxia*da2dzic
                  hessy(3,ic) = hessy(3,ic) + d2eda1a2*da1dyia*da2dzic
                  hessz(3,ic) = hessz(3,ic) + d2eda1a2*da1dzia*da2dzic
                  hessx(1,id) = hessx(1,id) + d2eda1a2*da1dxia*da2dxid
                  hessy(1,id) = hessy(1,id) + d2eda1a2*da1dyia*da2dxid
                  hessz(1,id) = hessz(1,id) + d2eda1a2*da1dzia*da2dxid
                  hessx(2,id) = hessx(2,id) + d2eda1a2*da1dxia*da2dyid
                  hessy(2,id) = hessy(2,id) + d2eda1a2*da1dyia*da2dyid
                  hessz(2,id) = hessz(2,id) + d2eda1a2*da1dzia*da2dyid
                  hessx(3,id) = hessx(3,id) + d2eda1a2*da1dxia*da2dzid
                  hessy(3,id) = hessy(3,id) + d2eda1a2*da1dyia*da2dzid
                  hessz(3,id) = hessz(3,id) + d2eda1a2*da1dzia*da2dzid
                  hessx(1,ie) = hessx(1,ie) + d2eda1a2*da1dxia*da2dxie
                  hessy(1,ie) = hessy(1,ie) + d2eda1a2*da1dyia*da2dxie
                  hessz(1,ie) = hessz(1,ie) + d2eda1a2*da1dzia*da2dxie
                  hessx(2,ie) = hessx(2,ie) + d2eda1a2*da1dxia*da2dyie
                  hessy(2,ie) = hessy(2,ie) + d2eda1a2*da1dyia*da2dyie
                  hessz(2,ie) = hessz(2,ie) + d2eda1a2*da1dzia*da2dyie
                  hessx(3,ie) = hessx(3,ie) + d2eda1a2*da1dxia*da2dzie
                  hessy(3,ie) = hessy(3,ie) + d2eda1a2*da1dyia*da2dzie
                  hessz(3,ie) = hessz(3,ie) + d2eda1a2*da1dzia*da2dzie
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedang2*dxibxib2
     &                             + d2eda2a2*da2dxib*da2dxib
     &                             + 2.0d0*d2eda1a2*da1dxib*da2dxib
                  hessy(1,ib) = hessy(1,ib) + dedang2*dxibyib2
     &                             + d2eda2a2*da2dxib*da2dyib
     &                             + d2eda1a2*da1dxib*da2dyib
     &                             + d2eda1a2*da2dxib*da1dyib
                  hessz(1,ib) = hessz(1,ib) + dedang2*dxibzib2
     &                             + d2eda2a2*da2dxib*da2dzib
     &                             + d2eda1a2*da1dxib*da2dzib
     &                             + d2eda1a2*da2dxib*da1dzib
                  hessx(2,ib) = hessx(2,ib) + dedang2*dxibyib2
     &                             + d2eda2a2*da2dxib*da2dyib
     &                             + d2eda1a2*da1dxib*da2dyib
     &                             + d2eda1a2*da2dxib*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang2*dyibyib2
     &                             + d2eda2a2*da2dyib*da2dyib
     &                             + d2eda1a2*da1dyib*da2dyib
     &                             + d2eda1a2*da2dyib*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang2*dyibzib2
     &                             + d2eda2a2*da2dyib*da2dzib
     &                             + d2eda1a2*da1dyib*da2dzib
     &                             + d2eda1a2*da2dyib*da1dzib
                  hessx(3,ib) = hessx(3,ib) + dedang2*dxibzib2
     &                             + d2eda2a2*da2dxib*da2dzib
     &                             + d2eda1a2*da1dxib*da2dzib
     &                             + d2eda1a2*da2dxib*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang2*dyibzib2
     &                             + d2eda2a2*da2dyib*da2dzib
     &                             + d2eda1a2*da1dyib*da2dzib
     &                             + d2eda1a2*da2dyib*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang2*dzibzib2
     &                             + d2eda2a2*da2dzib*da2dzib
     &                             + d2eda1a2*da1dzib*da2dzib
     &                             + d2eda1a2*da2dzib*da1dzib
                  hessx(1,ia) = hessx(1,ia) + d2eda1a2*da2dxib*da1dxia
                  hessy(1,ia) = hessy(1,ia) + d2eda1a2*da2dyib*da1dxia
                  hessz(1,ia) = hessz(1,ia) + d2eda1a2*da2dzib*da1dxia
                  hessx(2,ia) = hessx(2,ia) + d2eda1a2*da2dxib*da1dyia
                  hessy(2,ia) = hessy(2,ia) + d2eda1a2*da2dyib*da1dyia
                  hessz(2,ia) = hessz(2,ia) + d2eda1a2*da2dzib*da1dyia
                  hessx(3,ia) = hessx(3,ia) + d2eda1a2*da2dxib*da1dzia
                  hessy(3,ia) = hessy(3,ia) + d2eda1a2*da2dyib*da1dzia
                  hessz(3,ia) = hessz(3,ia) + d2eda1a2*da2dzib*da1dzia
                  hessx(1,ic) = hessx(1,ic) + dedang2*dxibxic2
     &                             + d2eda2a2*da2dxib*da2dxic
     &                             + d2eda1a2*da1dxib*da2dxic
     &                             + d2eda1a2*da2dxib*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang2*dyibxic2
     &                             + d2eda2a2*da2dyib*da2dxic
     &                             + d2eda1a2*da1dyib*da2dxic
     &                             + d2eda1a2*da2dyib*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang2*dzibxic2
     &                             + d2eda2a2*da2dzib*da2dxic
     &                             + d2eda1a2*da1dzib*da2dxic
     &                             + d2eda1a2*da2dzib*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang2*dxibyic2
     &                             + d2eda2a2*da2dxib*da2dyic
     &                             + d2eda1a2*da1dxib*da2dyic
     &                             + d2eda1a2*da2dxib*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang2*dyibyic2
     &                             + d2eda2a2*da2dyib*da2dyic
     &                             + d2eda1a2*da1dyib*da2dyic
     &                             + d2eda1a2*da2dyib*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang2*dzibyic2
     &                             + d2eda2a2*da2dzib*da2dyic
     &                             + d2eda1a2*da1dzib*da2dyic
     &                             + d2eda1a2*da2dzib*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang2*dxibzic2
     &                             + d2eda2a2*da2dxib*da2dzic
     &                             + d2eda1a2*da1dxib*da2dzic
     &                             + d2eda1a2*da2dxib*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang2*dyibzic2
     &                             + d2eda2a2*da2dyib*da2dzic
     &                             + d2eda1a2*da1dyib*da2dzic
     &                             + d2eda1a2*da2dyib*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang2*dzibzic2
     &                             + d2eda2a2*da2dzib*da2dzic
     &                             + d2eda1a2*da1dzib*da2dzic
     &                             + d2eda1a2*da2dzib*da1dzic
                  hessx(1,id) = hessx(1,id) + dedang2*dxibxid2
     &                             + d2eda2a2*da2dxib*da2dxid
     &                             + d2eda1a2*da1dxib*da2dxid
     &                             + d2eda1a2*da2dxib*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang2*dyibxid2
     &                             + d2eda2a2*da2dyib*da2dxid
     &                             + d2eda1a2*da1dyib*da2dxid
     &                             + d2eda1a2*da2dyib*da1dxid
                  hessz(1,id) = hessz(1,id) + dedang2*dzibxid2
     &                             + d2eda2a2*da2dzib*da2dxid
     &                             + d2eda1a2*da1dzib*da2dxid
     &                             + d2eda1a2*da2dzib*da1dxid
                  hessx(2,id) = hessx(2,id) + dedang2*dxibyid2
     &                             + d2eda2a2*da2dxib*da2dyid
     &                             + d2eda1a2*da1dxib*da2dyid
     &                             + d2eda1a2*da2dxib*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang2*dyibyid2
     &                             + d2eda2a2*da2dyib*da2dyid
     &                             + d2eda1a2*da1dyib*da2dyid
     &                             + d2eda1a2*da2dyib*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang2*dzibyid2
     &                             + d2eda2a2*da2dzib*da2dyid
     &                             + d2eda1a2*da1dzib*da2dyid
     &                             + d2eda1a2*da2dzib*da1dyid
                  hessx(3,id) = hessx(3,id) + dedang2*dxibzid2
     &                             + d2eda2a2*da2dxib*da2dzid
     &                             + d2eda1a2*da1dxib*da2dzid
     &                             + d2eda1a2*da2dxib*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang2*dyibzid2
     &                             + d2eda2a2*da2dyib*da2dzid
     &                             + d2eda1a2*da1dyib*da2dzid
     &                             + d2eda1a2*da2dyib*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang2*dzibzid2
     &                             + d2eda2a2*da2dzib*da2dzid
     &                             + d2eda1a2*da1dzib*da2dzid
     &                             + d2eda1a2*da2dzib*da1dzid
                  hessx(1,ie) = hessx(1,ie) + dedang2*dxibxie2
     &                             + d2eda2a2*da2dxib*da2dxie
     &                             + d2eda1a2*da1dxib*da2dxie
                  hessy(1,ie) = hessy(1,ie) + dedang2*dyibxie2
     &                             + d2eda2a2*da2dyib*da2dxie
     &                             + d2eda1a2*da1dyib*da2dxie
                  hessz(1,ie) = hessz(1,ie) + dedang2*dzibxie2
     &                             + d2eda2a2*da2dzib*da2dxie
     &                             + d2eda1a2*da1dzib*da2dxie
                  hessx(2,ie) = hessx(2,ie) + dedang2*dxibyie2
     &                             + d2eda2a2*da2dxib*da2dyie
     &                             + d2eda1a2*da1dxib*da2dyie
                  hessy(2,ie) = hessy(2,ie) + dedang2*dyibyie2
     &                             + d2eda2a2*da2dyib*da2dyie
     &                             + d2eda1a2*da1dyib*da2dyie
                  hessz(2,ie) = hessz(2,ie) + dedang2*dzibyie2
     &                             + d2eda2a2*da2dzib*da2dyie
     &                             + d2eda1a2*da1dzib*da2dyie
                  hessx(3,ie) = hessx(3,ie) + dedang2*dxibzie2
     &                             + d2eda2a2*da2dxib*da2dzie
     &                             + d2eda1a2*da1dxib*da2dzie
                  hessy(3,ie) = hessy(3,ie) + dedang2*dyibzie2
     &                             + d2eda2a2*da2dyib*da2dzie
     &                             + d2eda1a2*da1dyib*da2dzie
                  hessz(3,ie) = hessz(3,ie) + dedang2*dzibzie2
     &                             + d2eda2a2*da2dzib*da2dzie
     &                             + d2eda1a2*da1dzib*da2dzie
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedang2*dxicxic2
     &                             + d2eda2a2*da2dxic*da2dxic
     &                             + d2eda1a2*da1dxic*da2dxic
     &                             + d2eda1a2*da2dxic*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang2*dxicyic2
     &                             + d2eda2a2*da2dxic*da2dyic
     &                             + d2eda1a2*da1dxic*da2dyic
     &                             + d2eda1a2*da2dxic*da1dyic
                  hessz(1,ic) = hessz(1,ic) + dedang2*dxiczic2
     &                             + d2eda2a2*da2dxic*da2dzic
     &                             + d2eda1a2*da1dxic*da2dzic
     &                             + d2eda1a2*da2dxic*da1dzic
                  hessx(2,ic) = hessx(2,ic) + dedang2*dxicyic2
     &                             + d2eda2a2*da2dxic*da2dyic
     &                             + d2eda1a2*da1dxic*da2dyic
     &                             + d2eda1a2*da2dxic*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang2*dyicyic2
     &                             + d2eda2a2*da2dyic*da2dyic
     &                             + d2eda1a2*da1dyic*da2dyic
     &                             + d2eda1a2*da2dyic*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang2*dyiczic2
     &                             + d2eda2a2*da2dyic*da2dzic
     &                             + d2eda1a2*da1dyic*da2dzic
     &                             + d2eda1a2*da2dyic*da1dzic
                  hessx(3,ic) = hessx(3,ic) + dedang2*dxiczic2
     &                             + d2eda2a2*da2dxic*da2dzic
     &                             + d2eda1a2*da1dxic*da2dzic
     &                             + d2eda1a2*da2dxic*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang2*dyiczic2
     &                             + d2eda2a2*da2dyic*da2dzic
     &                             + d2eda1a2*da1dyic*da2dzic
     &                             + d2eda1a2*da2dyic*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang2*dziczic2
     &                             + d2eda2a2*da2dzic*da2dzic
     &                             + d2eda1a2*da1dzic*da2dzic
     &                             + d2eda1a2*da2dzic*da1dzic
                  hessx(1,ia) = hessx(1,ia) + d2eda1a2*da2dxic*da1dxia
                  hessy(1,ia) = hessy(1,ia) + d2eda1a2*da2dyic*da1dxia
                  hessz(1,ia) = hessz(1,ia) + d2eda1a2*da2dzic*da1dxia
                  hessx(2,ia) = hessx(2,ia) + d2eda1a2*da2dxic*da1dyia
                  hessy(2,ia) = hessy(2,ia) + d2eda1a2*da2dyic*da1dyia
                  hessz(2,ia) = hessz(2,ia) + d2eda1a2*da2dzic*da1dyia
                  hessx(3,ia) = hessx(3,ia) + d2eda1a2*da2dxic*da1dzia
                  hessy(3,ia) = hessy(3,ia) + d2eda1a2*da2dyic*da1dzia
                  hessz(3,ia) = hessz(3,ia) + d2eda1a2*da2dzic*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang2*dxibxic2
     &                             + d2eda2a2*da2dxic*da2dxib
     &                             + d2eda1a2*da1dxic*da2dxib
     &                             + d2eda1a2*da2dxic*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang2*dxibyic2
     &                             + d2eda2a2*da2dyic*da2dxib
     &                             + d2eda1a2*da1dyic*da2dxib
     &                             + d2eda1a2*da2dyic*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang2*dxibzic2
     &                             + d2eda2a2*da2dzic*da2dxib
     &                             + d2eda1a2*da1dzic*da2dxib
     &                             + d2eda1a2*da2dzic*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang2*dyibxic2
     &                             + d2eda2a2*da2dxic*da2dyib
     &                             + d2eda1a2*da1dxic*da2dyib
     &                             + d2eda1a2*da2dxic*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang2*dyibyic2
     &                             + d2eda2a2*da2dyic*da2dyib
     &                             + d2eda1a2*da1dyic*da2dyib
     &                             + d2eda1a2*da2dyic*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang2*dyibzic2
     &                             + d2eda2a2*da2dzic*da2dyib
     &                             + d2eda1a2*da1dzic*da2dyib
     &                             + d2eda1a2*da2dzic*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang2*dzibxic2
     &                             + d2eda2a2*da2dxic*da2dzib
     &                             + d2eda1a2*da1dxic*da2dzib
     &                             + d2eda1a2*da2dxic*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang2*dzibyic2
     &                             + d2eda2a2*da2dyic*da2dzib
     &                             + d2eda1a2*da1dyic*da2dzib
     &                             + d2eda1a2*da2dyic*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang2*dzibzic2
     &                             + d2eda2a2*da2dzic*da2dzib
     &                             + d2eda1a2*da1dzic*da2dzib
     &                             + d2eda1a2*da2dzic*da1dzib
                  hessx(1,id) = hessx(1,id) + dedang2*dxicxid2
     &                             + d2eda2a2*da2dxic*da2dxid
     &                             + d2eda1a2*da1dxic*da2dxid
     &                             + d2eda1a2*da2dxic*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang2*dyicxid2
     &                             + d2eda2a2*da2dyic*da2dxid
     &                             + d2eda1a2*da1dyic*da2dxid
     &                             + d2eda1a2*da2dyic*da1dxid
                  hessz(1,id) = hessz(1,id) + dedang2*dzicxid2
     &                             + d2eda2a2*da2dzic*da2dxid
     &                             + d2eda1a2*da1dzic*da2dxid
     &                             + d2eda1a2*da2dzic*da1dxid
                  hessx(2,id) = hessx(2,id) + dedang2*dxicyid2
     &                             + d2eda2a2*da2dxic*da2dyid
     &                             + d2eda1a2*da1dxic*da2dyid
     &                             + d2eda1a2*da2dxic*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang2*dyicyid2
     &                             + d2eda2a2*da2dyic*da2dyid
     &                             + d2eda1a2*da1dyic*da2dyid
     &                             + d2eda1a2*da2dyic*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang2*dzicyid2
     &                             + d2eda2a2*da2dzic*da2dyid
     &                             + d2eda1a2*da1dzic*da2dyid
     &                             + d2eda1a2*da2dzic*da1dyid
                  hessx(3,id) = hessx(3,id) + dedang2*dxiczid2
     &                             + d2eda2a2*da2dxic*da2dzid
     &                             + d2eda1a2*da1dxic*da2dzid
     &                             + d2eda1a2*da2dxic*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang2*dyiczid2
     &                             + d2eda2a2*da2dyic*da2dzid
     &                             + d2eda1a2*da1dyic*da2dzid
     &                             + d2eda1a2*da2dyic*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang2*dziczid2
     &                             + d2eda2a2*da2dzic*da2dzid
     &                             + d2eda1a2*da1dzic*da2dzid
     &                             + d2eda1a2*da2dzic*da1dzid
                  hessx(1,ie) = hessx(1,ie) + dedang2*dxicxie2
     &                             + d2eda2a2*da2dxic*da2dxie
     &                             + d2eda1a2*da1dxic*da2dxie
                  hessy(1,ie) = hessy(1,ie) + dedang2*dyicxie2
     &                             + d2eda2a2*da2dyic*da2dxie
     &                             + d2eda1a2*da1dyic*da2dxie
                  hessz(1,ie) = hessz(1,ie) + dedang2*dzicxie2
     &                             + d2eda2a2*da2dzic*da2dxie
     &                             + d2eda1a2*da1dzic*da2dxie
                  hessx(2,ie) = hessx(2,ie) + dedang2*dxicyie2
     &                             + d2eda2a2*da2dxic*da2dyie
     &                             + d2eda1a2*da1dxic*da2dyie
                  hessy(2,ie) = hessy(2,ie) + dedang2*dyicyie2
     &                             + d2eda2a2*da2dyic*da2dyie
     &                             + d2eda1a2*da1dyic*da2dyie
                  hessz(2,ie) = hessz(2,ie) + dedang2*dzicyie2
     &                             + d2eda2a2*da2dzic*da2dyie
     &                             + d2eda1a2*da1dzic*da2dyie
                  hessx(3,ie) = hessx(3,ie) + dedang2*dxiczie2
     &                             + d2eda2a2*da2dxic*da2dzie
     &                             + d2eda1a2*da1dxic*da2dzie
                  hessy(3,ie) = hessy(3,ie) + dedang2*dyiczie2
     &                             + d2eda2a2*da2dyic*da2dzie
     &                             + d2eda1a2*da1dyic*da2dzie
                  hessz(3,ie) = hessz(3,ie) + dedang2*dziczie2
     &                             + d2eda2a2*da2dzic*da2dzie
     &                             + d2eda1a2*da1dzic*da2dzie
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedang2*dxidxid2
     &                             + d2eda2a2*da2dxid*da2dxid
     &                             + d2eda1a2*da1dxid*da2dxid
     &                             + d2eda1a2*da2dxid*da1dxid
                  hessy(1,id) = hessy(1,id) + dedang2*dxidyid2
     &                             + d2eda2a2*da2dxid*da2dyid
     &                             + d2eda1a2*da1dxid*da2dyid
     &                             + d2eda1a2*da2dxid*da1dyid
                  hessz(1,id) = hessz(1,id) + dedang2*dxidzid2
     &                             + d2eda2a2*da2dxid*da2dzid
     &                             + d2eda1a2*da1dxid*da2dzid
     &                             + d2eda1a2*da2dxid*da1dzid
                  hessx(2,id) = hessx(2,id) + dedang2*dxidyid2
     &                             + d2eda2a2*da2dxid*da2dyid
     &                             + d2eda1a2*da1dxid*da2dyid
     &                             + d2eda1a2*da2dxid*da1dyid
                  hessy(2,id) = hessy(2,id) + dedang2*dyidyid2
     &                             + d2eda2a2*da2dyid*da2dyid
     &                             + d2eda1a2*da1dyid*da2dyid
     &                             + d2eda1a2*da2dyid*da1dyid
                  hessz(2,id) = hessz(2,id) + dedang2*dyidzid2
     &                             + d2eda2a2*da2dyid*da2dzid
     &                             + d2eda1a2*da1dyid*da2dzid
     &                             + d2eda1a2*da2dyid*da1dzid
                  hessx(3,id) = hessx(3,id) + dedang2*dxidzid2
     &                             + d2eda2a2*da2dxid*da2dzid
     &                             + d2eda1a2*da1dxid*da2dzid
     &                             + d2eda1a2*da2dxid*da1dzid
                  hessy(3,id) = hessy(3,id) + dedang2*dyidzid2
     &                             + d2eda2a2*da2dyid*da2dzid
     &                             + d2eda1a2*da1dyid*da2dzid
     &                             + d2eda1a2*da2dyid*da1dzid
                  hessz(3,id) = hessz(3,id) + dedang2*dzidzid2
     &                             + d2eda2a2*da2dzid*da2dzid
     &                             + d2eda1a2*da1dzid*da2dzid
     &                             + d2eda1a2*da2dzid*da1dzid
                  hessx(1,ia) = hessx(1,ia) + d2eda1a2*da2dxid*da1dxia
                  hessy(1,ia) = hessy(1,ia) + d2eda1a2*da2dyid*da1dxia
                  hessz(1,ia) = hessz(1,ia) + d2eda1a2*da2dzid*da1dxia
                  hessx(2,ia) = hessx(2,ia) + d2eda1a2*da2dxid*da1dyia
                  hessy(2,ia) = hessy(2,ia) + d2eda1a2*da2dyid*da1dyia
                  hessz(2,ia) = hessz(2,ia) + d2eda1a2*da2dzid*da1dyia
                  hessx(3,ia) = hessx(3,ia) + d2eda1a2*da2dxid*da1dzia
                  hessy(3,ia) = hessy(3,ia) + d2eda1a2*da2dyid*da1dzia
                  hessz(3,ia) = hessz(3,ia) + d2eda1a2*da2dzid*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang2*dxibxid2
     &                             + d2eda2a2*da2dxid*da2dxib
     &                             + d2eda1a2*da1dxid*da2dxib
     &                             + d2eda1a2*da2dxid*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang2*dxibyid2
     &                             + d2eda2a2*da2dyid*da2dxib
     &                             + d2eda1a2*da1dyid*da2dxib
     &                             + d2eda1a2*da2dyid*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang2*dxibzid2
     &                             + d2eda2a2*da2dzid*da2dxib
     &                             + d2eda1a2*da1dzid*da2dxib
     &                             + d2eda1a2*da2dzid*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang2*dyibxid2
     &                             + d2eda2a2*da2dxid*da2dyib
     &                             + d2eda1a2*da1dxid*da2dyib
     &                             + d2eda1a2*da2dxid*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang2*dyibyid2
     &                             + d2eda2a2*da2dyid*da2dyib
     &                             + d2eda1a2*da1dyid*da2dyib
     &                             + d2eda1a2*da2dyid*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang2*dyibzid2
     &                             + d2eda2a2*da2dzid*da2dyib
     &                             + d2eda1a2*da1dzid*da2dyib
     &                             + d2eda1a2*da2dzid*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang2*dzibxid2
     &                             + d2eda2a2*da2dxid*da2dzib
     &                             + d2eda1a2*da1dxid*da2dzib
     &                             + d2eda1a2*da2dxid*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang2*dzibyid2
     &                             + d2eda2a2*da2dyid*da2dzib
     &                             + d2eda1a2*da1dyid*da2dzib
     &                             + d2eda1a2*da2dyid*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang2*dzibzid2
     &                             + d2eda2a2*da2dzid*da2dzib
     &                             + d2eda1a2*da1dzid*da2dzib
     &                             + d2eda1a2*da2dzid*da1dzib
                  hessx(1,ic) = hessx(1,ic) + dedang2*dxicxid2
     &                             + d2eda2a2*da2dxid*da2dxic
     &                             + d2eda1a2*da1dxid*da2dxic
     &                             + d2eda1a2*da2dxid*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang2*dxicyid2
     &                             + d2eda2a2*da2dyid*da2dxic
     &                             + d2eda1a2*da1dyid*da2dxic
     &                             + d2eda1a2*da2dyid*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang2*dxiczid2
     &                             + d2eda2a2*da2dzid*da2dxic
     &                             + d2eda1a2*da1dzid*da2dxic
     &                             + d2eda1a2*da2dzid*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang2*dyicxid2
     &                             + d2eda2a2*da2dxid*da2dyic
     &                             + d2eda1a2*da1dxid*da2dyic
     &                             + d2eda1a2*da2dxid*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang2*dyicyid2
     &                             + d2eda2a2*da2dyid*da2dyic
     &                             + d2eda1a2*da1dyid*da2dyic
     &                             + d2eda1a2*da2dyid*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang2*dyiczid2
     &                             + d2eda2a2*da2dzid*da2dyic
     &                             + d2eda1a2*da1dzid*da2dyic
     &                             + d2eda1a2*da2dzid*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang2*dzicxid2
     &                             + d2eda2a2*da2dxid*da2dzic
     &                             + d2eda1a2*da1dxid*da2dzic
     &                             + d2eda1a2*da2dxid*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang2*dzicyid2
     &                             + d2eda2a2*da2dyid*da2dzic
     &                             + d2eda1a2*da1dyid*da2dzic
     &                             + d2eda1a2*da2dyid*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang2*dziczid2
     &                             + d2eda2a2*da2dzid*da2dzic
     &                             + d2eda1a2*da1dzid*da2dzic
     &                             + d2eda1a2*da2dzid*da1dzic
                  hessx(1,ie) = hessx(1,ie) + dedang2*dxidxie2
     &                             + d2eda2a2*da2dxid*da2dxie
     &                             + d2eda1a2*da1dxid*da2dxie
                  hessy(1,ie) = hessy(1,ie) + dedang2*dyidxie2
     &                             + d2eda2a2*da2dyid*da2dxie
     &                             + d2eda1a2*da1dyid*da2dxie
                  hessz(1,ie) = hessz(1,ie) + dedang2*dzidxie2
     &                             + d2eda2a2*da2dzid*da2dxie
     &                             + d2eda1a2*da1dzid*da2dxie
                  hessx(2,ie) = hessx(2,ie) + dedang2*dxidyie2
     &                             + d2eda2a2*da2dxid*da2dyie
     &                             + d2eda1a2*da1dxid*da2dyie
                  hessy(2,ie) = hessy(2,ie) + dedang2*dyidyie2
     &                             + d2eda2a2*da2dyid*da2dyie
     &                             + d2eda1a2*da1dyid*da2dyie
                  hessz(2,ie) = hessz(2,ie) + dedang2*dzidyie2
     &                             + d2eda2a2*da2dzid*da2dyie
     &                             + d2eda1a2*da1dzid*da2dyie
                  hessx(3,ie) = hessx(3,ie) + dedang2*dxidzie2
     &                             + d2eda2a2*da2dxid*da2dzie
     &                             + d2eda1a2*da1dxid*da2dzie
                  hessy(3,ie) = hessy(3,ie) + dedang2*dyidzie2
     &                             + d2eda2a2*da2dyid*da2dzie
     &                             + d2eda1a2*da1dyid*da2dzie
                  hessz(3,ie) = hessz(3,ie) + dedang2*dzidzie2
     &                             + d2eda2a2*da2dzid*da2dzie
     &                             + d2eda1a2*da1dzid*da2dzie
               else if (i .eq. ie) then
                  hessx(1,ie) = hessx(1,ie) + dedang2*dxiexie2
     &                             + d2eda2a2*da2dxie*da2dxie
                  hessy(1,ie) = hessy(1,ie) + dedang2*dxieyie2
     &                             + d2eda2a2*da2dxie*da2dyie
                  hessz(1,ie) = hessz(1,ie) + dedang2*dxiezie2
     &                             + d2eda2a2*da2dxie*da2dzie
                  hessx(2,ie) = hessx(2,ie) + dedang2*dxieyie2
     &                             + d2eda2a2*da2dxie*da2dyie
                  hessy(2,ie) = hessy(2,ie) + dedang2*dyieyie2
     &                             + d2eda2a2*da2dyie*da2dyie
                  hessz(2,ie) = hessz(2,ie) + dedang2*dyiezie2
     &                             + d2eda2a2*da2dyie*da2dzie
                  hessx(3,ie) = hessx(3,ie) + dedang2*dxiezie2
     &                             + d2eda2a2*da2dxie*da2dzie
                  hessy(3,ie) = hessy(3,ie) + dedang2*dyiezie2
     &                             + d2eda2a2*da2dyie*da2dzie
                  hessz(3,ie) = hessz(3,ie) + dedang2*dziezie2
     &                             + d2eda2a2*da2dzie*da2dzie
                  hessx(1,ia) = hessx(1,ia) + d2eda1a2*da2dxie*da1dxia
                  hessy(1,ia) = hessy(1,ia) + d2eda1a2*da2dyie*da1dxia
                  hessz(1,ia) = hessz(1,ia) + d2eda1a2*da2dzie*da1dxia
                  hessx(2,ia) = hessx(2,ia) + d2eda1a2*da2dxie*da1dyia
                  hessy(2,ia) = hessy(2,ia) + d2eda1a2*da2dyie*da1dyia
                  hessz(2,ia) = hessz(2,ia) + d2eda1a2*da2dzie*da1dyia
                  hessx(3,ia) = hessx(3,ia) + d2eda1a2*da2dxie*da1dzia
                  hessy(3,ia) = hessy(3,ia) + d2eda1a2*da2dyie*da1dzia
                  hessz(3,ia) = hessz(3,ia) + d2eda1a2*da2dzie*da1dzia
                  hessx(1,ib) = hessx(1,ib) + dedang2*dxibxie2
     &                             + d2eda2a2*da2dxie*da2dxib
     &                             + d2eda1a2*da2dxie*da1dxib
                  hessy(1,ib) = hessy(1,ib) + dedang2*dxibyie2
     &                             + d2eda2a2*da2dyie*da2dxib
     &                             + d2eda1a2*da2dyie*da1dxib
                  hessz(1,ib) = hessz(1,ib) + dedang2*dxibzie2
     &                             + d2eda2a2*da2dzie*da2dxib
     &                             + d2eda1a2*da2dzie*da1dxib
                  hessx(2,ib) = hessx(2,ib) + dedang2*dyibxie2
     &                             + d2eda2a2*da2dxie*da2dyib
     &                             + d2eda1a2*da2dxie*da1dyib
                  hessy(2,ib) = hessy(2,ib) + dedang2*dyibyie2
     &                             + d2eda2a2*da2dyie*da2dyib
     &                             + d2eda1a2*da2dyie*da1dyib
                  hessz(2,ib) = hessz(2,ib) + dedang2*dyibzie2
     &                             + d2eda2a2*da2dzie*da2dyib
     &                             + d2eda1a2*da2dzie*da1dyib
                  hessx(3,ib) = hessx(3,ib) + dedang2*dzibxie2
     &                             + d2eda2a2*da2dxie*da2dzib
     &                             + d2eda1a2*da2dxie*da1dzib
                  hessy(3,ib) = hessy(3,ib) + dedang2*dzibyie2
     &                             + d2eda2a2*da2dyie*da2dzib
     &                             + d2eda1a2*da2dyie*da1dzib
                  hessz(3,ib) = hessz(3,ib) + dedang2*dzibzie2
     &                             + d2eda2a2*da2dzie*da2dzib
     &                             + d2eda1a2*da2dzie*da1dzib
                  hessx(1,ic) = hessx(1,ic) + dedang2*dxicxie2
     &                             + d2eda2a2*da2dxie*da2dxic
     &                             + d2eda1a2*da2dxie*da1dxic
                  hessy(1,ic) = hessy(1,ic) + dedang2*dxicyie2
     &                             + d2eda2a2*da2dyie*da2dxic
     &                             + d2eda1a2*da2dyie*da1dxic
                  hessz(1,ic) = hessz(1,ic) + dedang2*dxiczie2
     &                             + d2eda2a2*da2dzie*da2dxic
     &                             + d2eda1a2*da2dzie*da1dxic
                  hessx(2,ic) = hessx(2,ic) + dedang2*dyicxie2
     &                             + d2eda2a2*da2dxie*da2dyic
     &                             + d2eda1a2*da2dxie*da1dyic
                  hessy(2,ic) = hessy(2,ic) + dedang2*dyicyie2
     &                             + d2eda2a2*da2dyie*da2dyic
     &                             + d2eda1a2*da2dyie*da1dyic
                  hessz(2,ic) = hessz(2,ic) + dedang2*dyiczie2
     &                             + d2eda2a2*da2dzie*da2dyic
     &                             + d2eda1a2*da2dzie*da1dyic
                  hessx(3,ic) = hessx(3,ic) + dedang2*dzicxie2
     &                             + d2eda2a2*da2dxie*da2dzic
     &                             + d2eda1a2*da2dxie*da1dzic
                  hessy(3,ic) = hessy(3,ic) + dedang2*dzicyie2
     &                             + d2eda2a2*da2dyie*da2dzic
     &                             + d2eda1a2*da2dyie*da1dzic
                  hessz(3,ic) = hessz(3,ic) + dedang2*dziczie2
     &                             + d2eda2a2*da2dzie*da2dzic
     &                             + d2eda1a2*da2dzie*da1dzic
                  hessx(1,id) = hessx(1,id) + dedang2*dxidxie2
     &                             + d2eda2a2*da2dxid*da2dxie
     &                             + d2eda1a2*da1dxid*da2dxie
                  hessy(1,id) = hessy(1,id) + dedang2*dxidyie2
     &                             + d2eda2a2*da2dxid*da2dyie
     &                             + d2eda1a2*da1dxid*da2dyie
                  hessz(1,id) = hessz(1,id) + dedang2*dxidzie2
     &                             + d2eda2a2*da2dxid*da2dzie
     &                             + d2eda1a2*da1dxid*da2dzie
                  hessx(2,id) = hessx(2,id) + dedang2*dyidxie2
     &                             + d2eda2a2*da2dyid*da2dxie
     &                             + d2eda1a2*da1dyid*da2dxie
                  hessy(2,id) = hessy(2,id) + dedang2*dyidyie2
     &                             + d2eda2a2*da2dyid*da2dyie
     &                             + d2eda1a2*da1dyid*da2dyie
                  hessz(2,id) = hessz(2,id) + dedang2*dyidzie2
     &                             + d2eda2a2*da2dyid*da2dzie
     &                             + d2eda1a2*da1dyid*da2dzie
                  hessx(3,id) = hessx(3,id) + dedang2*dzidxie2
     &                             + d2eda2a2*da2dzid*da2dxie
     &                             + d2eda1a2*da1dzid*da2dxie
                  hessy(3,id) = hessy(3,id) + dedang2*dzidyie2
     &                             + d2eda2a2*da2dzid*da2dyie
     &                             + d2eda1a2*da1dzid*da2dyie
                  hessz(3,id) = hessz(3,id) + dedang2*dzidzie2
     &                             + d2eda2a2*da2dzid*da2dzie
     &                             + d2eda1a2*da1dzid*da2dzie
               end if
            end if
         end if
      end do
      return
      end
