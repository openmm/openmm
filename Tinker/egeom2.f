c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine egeom2  --  atom-by-atom restraint Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "egeom2" calculates second derivatives of restraints
c     on positions, distances, angles and torsions as well
c     as Gaussian basin and spherical droplet restraints
c
c     note that the Hessian is discontinuous when an upper and
c     lower bound range is used instead of a single distance
c
c
      subroutine egeom2 (i)
      use sizes
      use atomid
      use atoms
      use bound
      use deriv
      use group
      use hessn
      use math
      use molcul
      use restrn
      implicit none
      integer i,j,k,m
      integer ia,ib,ic,id
      integer kpos,kdist,kang
      integer ktors,kchir
      real*8 xr,yr,zr,fgrp
      real*8 target,force
      real*8 dot,angle
      real*8 cosine,sine
      real*8 dt,dt2,deddt
      real*8 term,terma,termc
      real*8 termx,termy,termz
      real*8 de,d2eddt2
      real*8 d2e(3,3)
      real*8 dedphi,d2edphi2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xrab,yrab,zrab
      real*8 xrcb,yrcb,zrcb
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 rab2,rcb2
      real*8 xpo,ypo,zpo
      real*8 xp,yp,zp,rp,rp2
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 df1,df2,af1,af2
      real*8 tf1,tf2,t1,t2
      real*8 gf1,gf2
      real*8 weigh,ratio
      real*8 weigha,weighb
      real*8 xcm,ycm,zcm
      real*8 cf1,cf2,vol
      real*8 c1,c2,c3
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
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
      real*8 dxibxia,dxibyia,dxibzia
      real*8 dyibxia,dyibyia,dyibzia
      real*8 dzibxia,dzibyia,dzibzia
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid
      real*8 dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid
      real*8 dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid
      real*8 dzicxid,dzicyid,dziczid
      real*8 ddtdxid,ddtdyid,ddtdzid
      real*8 dedr,d2edr2,expterm
      real*8 xi,yi,zi,ri,ri2
      real*8 r,r2,r6,r12
      real*8 a,b,buffer
      logical proceed,intermol,linear
c
c
c     compute the Hessian elements for position restraints
c
      do kpos = 1, npfix
         ia = ipfix(kpos)
         proceed = (i .eq. ia)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,0,0,0,0,0)
         if (proceed) then
            xr = 0.0d0
            yr = 0.0d0
            zr = 0.0d0
            if (kpfix(1,i) .ne. 0)  xr = x(ia) - xpfix(kpos)
            if (kpfix(2,i) .ne. 0)  yr = y(ia) - ypfix(kpos)
            if (kpfix(3,i) .ne. 0)  zr = z(ia) - zpfix(kpos)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            force = pfix(1,kpos)
            dt = max(0.0d0,r-pfix(2,kpos))
            dt2 = dt * dt
            deddt = 2.0d0 * force
c
c     scale the interaction based on its group membership
c
            if (use_group)  deddt = deddt * fgrp
c
c     set the chain rule terms for the Hessian elements
c
            if (r .eq. 0.0d0) then
               de = deddt
               term = 0.0d0
            else
               de = deddt * dt/r
               term = (deddt-de) / r2
            end if
            termx = term * xr
            termy = term * yr
            termz = term * zr
            d2e(1,1) = termx*xr + de
            d2e(1,2) = termx*yr
            d2e(1,3) = termx*zr
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yr + de
            d2e(2,3) = termy*zr
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
            end do
         end if
      end do
c
c     compute the Hessian elements for distance restraints
c
      do kdist = 1, ndfix
         ia = idfix(1,kdist)
         ib = idfix(2,kdist)
         proceed = (i.eq.ia .or. i.eq.ib)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,0,0,0,0)
         if (proceed) then
            if (i .eq. ib) then
               ib = ia
               ia = i
            end if
            xr = x(ia) - x(ib)
            yr = y(ia) - y(ib)
            zr = z(ia) - z(ib)
            intermol = (molcule(ia) .ne. molcule(ib))
            if (use_bounds .and. intermol)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            force = dfix(1,kdist)
            df1 = dfix(2,kdist)
            df2 = dfix(3,kdist)
            target = r
            if (r .lt. df1)  target = df1
            if (r .gt. df2)  target = df2
            dt = r - target
            deddt = 2.0d0 * force
c
c     scale the interaction based on its group membership
c
            if (use_group)  deddt = deddt * fgrp
c
c     set the chain rule terms for the Hessian elements
c
            if (r .eq. 0.0d0) then
               r = 0.0001d0
               r2 = r * r
            end if
            de = deddt * dt/r
            term = (deddt-de) / r2
            termx = term * xr
            termy = term * yr
            termz = term * zr
            d2e(1,1) = termx*xr + de
            d2e(1,2) = termx*yr
            d2e(1,3) = termx*zr
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yr + de
            d2e(2,3) = termy*zr
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ib) = hessx(j,ib) - d2e(1,j)
               hessy(j,ib) = hessy(j,ib) - d2e(2,j)
               hessz(j,ib) = hessz(j,ib) - d2e(3,j)
            end do
         end if
      end do
c
c     compute the Hessian elements for angle restraints
c
      do kang = 1, nafix
         ia = iafix(1,kang)
         ib = iafix(2,kang)
         ic = iafix(3,kang)
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,0,0,0)
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
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               force = afix(1,kang)
               af1 = afix(2,kang)
               af2 = afix(3,kang)
               target = angle
               if (angle .lt. af1)  target = af1
               if (angle .gt. af2)  target = af2
               dt = angle - target
               dt2 = dt * dt
               deddt = 2.0d0 * force * dt * radian
               d2eddt2 = 2.0d0 * force * radian**2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  deddt = deddt * fgrp
                  d2eddt2 = d2eddt2 * fgrp
               end if
c
c     construct an orthogonal direction for linear angles
c
               linear = .false.
               if (rp .lt. 0.000001d0) then
                  linear = .true.
                  if (xab.ne.0.0d0 .and. yab.ne.0.0d0) then
                     xp = -yab
                     yp = xab
                     zp = 0.0d0
                  else if (xab.eq.0.0d0 .and. yab.eq.0.0d0) then
                     xp = 1.0d0
                     yp = 0.0d0
                     zp = 0.0d0
                  else if (xab.ne.0.0d0 .and. yab.eq.0.0d0) then
                     xp = 0.0d0
                     yp = 1.0d0
                     zp = 0.0d0
                  else if (xab.eq.0.0d0 .and. yab.ne.0.0d0) then
                     xp = 1.0d0
                     yp = 0.0d0
                     zp = 0.0d0
                  end if
                  rp = sqrt(xp*xp + yp*yp + zp*zp)
               end if
c
c     first derivatives of bond angle with respect to coordinates
c
   10          continue
               terma = -1.0d0 / (rab2*rp)
               termc = 1.0d0 / (rcb2*rp)
               ddtdxia = terma * (yab*zp-zab*yp)
               ddtdyia = terma * (zab*xp-xab*zp)
               ddtdzia = terma * (xab*yp-yab*xp)
               ddtdxic = termc * (ycb*zp-zcb*yp)
               ddtdyic = termc * (zcb*xp-xcb*zp)
               ddtdzic = termc * (xcb*yp-ycb*xp)
               ddtdxib = -ddtdxia - ddtdxic
               ddtdyib = -ddtdyia - ddtdyic
               ddtdzib = -ddtdzia - ddtdzic
c
c     abbreviations used in defining chain rule terms
c
               xrab = 2.0d0 * xab / rab2
               yrab = 2.0d0 * yab / rab2
               zrab = 2.0d0 * zab / rab2
               xrcb = 2.0d0 * xcb / rcb2
               yrcb = 2.0d0 * ycb / rcb2
               zrcb = 2.0d0 * zcb / rcb2
               rp2 = 1.0d0 / (rp*rp)
               xabp = (yab*zp-zab*yp) * rp2
               yabp = (zab*xp-xab*zp) * rp2
               zabp = (xab*yp-yab*xp) * rp2
               xcbp = (ycb*zp-zcb*yp) * rp2
               ycbp = (zcb*xp-xcb*zp) * rp2
               zcbp = (xcb*yp-ycb*xp) * rp2
c
c     chain rule terms for second derivative components
c
               dxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab)
               dxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab)
               dxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab)
               dyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab)
               dyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab)
               dziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab)
               dxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb)
               dxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb)
               dxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb)
               dyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb)
               dyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb)
               dziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb)
               dxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp
               dxiayic = -terma*xab*yab - ddtdxia*yabp
               dxiazic = -terma*xab*zab - ddtdxia*zabp
               dyiaxic = -terma*xab*yab - ddtdyia*xabp
               dyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp
               dyiazic = -terma*yab*zab - ddtdyia*zabp
               dziaxic = -terma*xab*zab - ddtdzia*xabp
               dziayic = -terma*yab*zab - ddtdzia*yabp
               dziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp
c
c     get some second derivative chain rule terms by difference
c
               dxibxia = -dxiaxia - dxiaxic
               dxibyia = -dxiayia - dyiaxic
               dxibzia = -dxiazia - dziaxic
               dyibxia = -dxiayia - dxiayic
               dyibyia = -dyiayia - dyiayic
               dyibzia = -dyiazia - dziayic
               dzibxia = -dxiazia - dxiazic
               dzibyia = -dyiazia - dyiazic
               dzibzia = -dziazia - dziazic
               dxibxic = -dxicxic - dxiaxic
               dxibyic = -dxicyic - dxiayic
               dxibzic = -dxiczic - dxiazic
               dyibxic = -dxicyic - dyiaxic
               dyibyic = -dyicyic - dyiayic
               dyibzic = -dyiczic - dyiazic
               dzibxic = -dxiczic - dziaxic
               dzibyic = -dyiczic - dziayic
               dzibzic = -dziczic - dziazic
               dxibxib = -dxibxia - dxibxic
               dxibyib = -dxibyia - dxibyic
               dxibzib = -dxibzia - dxibzic
               dyibyib = -dyibyia - dyibyic
               dyibzib = -dyibzia - dyibzic
               dzibzib = -dzibzia - dzibzic
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (ia .eq. i) then
                  hessx(1,ia) = hessx(1,ia) + deddt*dxiaxia
     &                               + d2eddt2*ddtdxia*ddtdxia
                  hessx(2,ia) = hessx(2,ia) + deddt*dxiayia
     &                               + d2eddt2*ddtdxia*ddtdyia
                  hessx(3,ia) = hessx(3,ia) + deddt*dxiazia
     &                               + d2eddt2*ddtdxia*ddtdzia
                  hessy(1,ia) = hessy(1,ia) + deddt*dxiayia
     &                               + d2eddt2*ddtdyia*ddtdxia
                  hessy(2,ia) = hessy(2,ia) + deddt*dyiayia
     &                               + d2eddt2*ddtdyia*ddtdyia
                  hessy(3,ia) = hessy(3,ia) + deddt*dyiazia
     &                               + d2eddt2*ddtdyia*ddtdzia
                  hessz(1,ia) = hessz(1,ia) + deddt*dxiazia
     &                               + d2eddt2*ddtdzia*ddtdxia
                  hessz(2,ia) = hessz(2,ia) + deddt*dyiazia
     &                               + d2eddt2*ddtdzia*ddtdyia
                  hessz(3,ia) = hessz(3,ia) + deddt*dziazia
     &                               + d2eddt2*ddtdzia*ddtdzia
                  hessx(1,ib) = hessx(1,ib) + deddt*dxibxia
     &                               + d2eddt2*ddtdxia*ddtdxib
                  hessx(2,ib) = hessx(2,ib) + deddt*dyibxia
     &                               + d2eddt2*ddtdxia*ddtdyib
                  hessx(3,ib) = hessx(3,ib) + deddt*dzibxia
     &                               + d2eddt2*ddtdxia*ddtdzib
                  hessy(1,ib) = hessy(1,ib) + deddt*dxibyia
     &                               + d2eddt2*ddtdyia*ddtdxib
                  hessy(2,ib) = hessy(2,ib) + deddt*dyibyia
     &                               + d2eddt2*ddtdyia*ddtdyib
                  hessy(3,ib) = hessy(3,ib) + deddt*dzibyia
     &                               + d2eddt2*ddtdyia*ddtdzib
                  hessz(1,ib) = hessz(1,ib) + deddt*dxibzia
     &                               + d2eddt2*ddtdzia*ddtdxib
                  hessz(2,ib) = hessz(2,ib) + deddt*dyibzia
     &                               + d2eddt2*ddtdzia*ddtdyib
                  hessz(3,ib) = hessz(3,ib) + deddt*dzibzia
     &                               + d2eddt2*ddtdzia*ddtdzib
                  hessx(1,ic) = hessx(1,ic) + deddt*dxiaxic
     &                               + d2eddt2*ddtdxia*ddtdxic
                  hessx(2,ic) = hessx(2,ic) + deddt*dxiayic
     &                               + d2eddt2*ddtdxia*ddtdyic
                  hessx(3,ic) = hessx(3,ic) + deddt*dxiazic
     &                               + d2eddt2*ddtdxia*ddtdzic
                  hessy(1,ic) = hessy(1,ic) + deddt*dyiaxic
     &                               + d2eddt2*ddtdyia*ddtdxic
                  hessy(2,ic) = hessy(2,ic) + deddt*dyiayic
     &                               + d2eddt2*ddtdyia*ddtdyic
                  hessy(3,ic) = hessy(3,ic) + deddt*dyiazic
     &                               + d2eddt2*ddtdyia*ddtdzic
                  hessz(1,ic) = hessz(1,ic) + deddt*dziaxic
     &                               + d2eddt2*ddtdzia*ddtdxic
                  hessz(2,ic) = hessz(2,ic) + deddt*dziayic
     &                               + d2eddt2*ddtdzia*ddtdyic
                  hessz(3,ic) = hessz(3,ic) + deddt*dziazic
     &                               + d2eddt2*ddtdzia*ddtdzic
               else if (ib .eq. i) then
                  hessx(1,ib) = hessx(1,ib) + deddt*dxibxib
     &                               + d2eddt2*ddtdxib*ddtdxib
                  hessx(2,ib) = hessx(2,ib) + deddt*dxibyib
     &                               + d2eddt2*ddtdxib*ddtdyib
                  hessx(3,ib) = hessx(3,ib) + deddt*dxibzib
     &                               + d2eddt2*ddtdxib*ddtdzib
                  hessy(1,ib) = hessy(1,ib) + deddt*dxibyib
     &                               + d2eddt2*ddtdyib*ddtdxib
                  hessy(2,ib) = hessy(2,ib) + deddt*dyibyib
     &                               + d2eddt2*ddtdyib*ddtdyib
                  hessy(3,ib) = hessy(3,ib) + deddt*dyibzib
     &                               + d2eddt2*ddtdyib*ddtdzib
                  hessz(1,ib) = hessz(1,ib) + deddt*dxibzib
     &                               + d2eddt2*ddtdzib*ddtdxib
                  hessz(2,ib) = hessz(2,ib) + deddt*dyibzib
     &                               + d2eddt2*ddtdzib*ddtdyib
                  hessz(3,ib) = hessz(3,ib) + deddt*dzibzib
     &                               + d2eddt2*ddtdzib*ddtdzib
                  hessx(1,ia) = hessx(1,ia) + deddt*dxibxia
     &                               + d2eddt2*ddtdxib*ddtdxia
                  hessx(2,ia) = hessx(2,ia) + deddt*dxibyia
     &                               + d2eddt2*ddtdxib*ddtdyia
                  hessx(3,ia) = hessx(3,ia) + deddt*dxibzia
     &                               + d2eddt2*ddtdxib*ddtdzia
                  hessy(1,ia) = hessy(1,ia) + deddt*dyibxia
     &                               + d2eddt2*ddtdyib*ddtdxia
                  hessy(2,ia) = hessy(2,ia) + deddt*dyibyia
     &                               + d2eddt2*ddtdyib*ddtdyia
                  hessy(3,ia) = hessy(3,ia) + deddt*dyibzia
     &                               + d2eddt2*ddtdyib*ddtdzia
                  hessz(1,ia) = hessz(1,ia) + deddt*dzibxia
     &                               + d2eddt2*ddtdzib*ddtdxia
                  hessz(2,ia) = hessz(2,ia) + deddt*dzibyia
     &                               + d2eddt2*ddtdzib*ddtdyia
                  hessz(3,ia) = hessz(3,ia) + deddt*dzibzia
     &                               + d2eddt2*ddtdzib*ddtdzia
                  hessx(1,ic) = hessx(1,ic) + deddt*dxibxic
     &                               + d2eddt2*ddtdxib*ddtdxic
                  hessx(2,ic) = hessx(2,ic) + deddt*dxibyic
     &                               + d2eddt2*ddtdxib*ddtdyic
                  hessx(3,ic) = hessx(3,ic) + deddt*dxibzic
     &                               + d2eddt2*ddtdxib*ddtdzic
                  hessy(1,ic) = hessy(1,ic) + deddt*dyibxic
     &                               + d2eddt2*ddtdyib*ddtdxic
                  hessy(2,ic) = hessy(2,ic) + deddt*dyibyic
     &                               + d2eddt2*ddtdyib*ddtdyic
                  hessy(3,ic) = hessy(3,ic) + deddt*dyibzic
     &                               + d2eddt2*ddtdyib*ddtdzic
                  hessz(1,ic) = hessz(1,ic) + deddt*dzibxic
     &                               + d2eddt2*ddtdzib*ddtdxic
                  hessz(2,ic) = hessz(2,ic) + deddt*dzibyic
     &                               + d2eddt2*ddtdzib*ddtdyic
                  hessz(3,ic) = hessz(3,ic) + deddt*dzibzic
     &                               + d2eddt2*ddtdzib*ddtdzic
               else if (ic .eq. i) then
                  hessx(1,ic) = hessx(1,ic) + deddt*dxicxic
     &                               + d2eddt2*ddtdxic*ddtdxic
                  hessx(2,ic) = hessx(2,ic) + deddt*dxicyic
     &                               + d2eddt2*ddtdxic*ddtdyic
                  hessx(3,ic) = hessx(3,ic) + deddt*dxiczic
     &                               + d2eddt2*ddtdxic*ddtdzic
                  hessy(1,ic) = hessy(1,ic) + deddt*dxicyic
     &                               + d2eddt2*ddtdyic*ddtdxic
                  hessy(2,ic) = hessy(2,ic) + deddt*dyicyic
     &                               + d2eddt2*ddtdyic*ddtdyic
                  hessy(3,ic) = hessy(3,ic) + deddt*dyiczic
     &                               + d2eddt2*ddtdyic*ddtdzic
                  hessz(1,ic) = hessz(1,ic) + deddt*dxiczic
     &                               + d2eddt2*ddtdzic*ddtdxic
                  hessz(2,ic) = hessz(2,ic) + deddt*dyiczic
     &                               + d2eddt2*ddtdzic*ddtdyic
                  hessz(3,ic) = hessz(3,ic) + deddt*dziczic
     &                               + d2eddt2*ddtdzic*ddtdzic
                  hessx(1,ib) = hessx(1,ib) + deddt*dxibxic
     &                               + d2eddt2*ddtdxic*ddtdxib
                  hessx(2,ib) = hessx(2,ib) + deddt*dyibxic
     &                               + d2eddt2*ddtdxic*ddtdyib
                  hessx(3,ib) = hessx(3,ib) + deddt*dzibxic
     &                               + d2eddt2*ddtdxic*ddtdzib
                  hessy(1,ib) = hessy(1,ib) + deddt*dxibyic
     &                               + d2eddt2*ddtdyic*ddtdxib
                  hessy(2,ib) = hessy(2,ib) + deddt*dyibyic
     &                               + d2eddt2*ddtdyic*ddtdyib
                  hessy(3,ib) = hessy(3,ib) + deddt*dzibyic
     &                               + d2eddt2*ddtdyic*ddtdzib
                  hessz(1,ib) = hessz(1,ib) + deddt*dxibzic
     &                               + d2eddt2*ddtdzic*ddtdxib
                  hessz(2,ib) = hessz(2,ib) + deddt*dyibzic
     &                               + d2eddt2*ddtdzic*ddtdyib
                  hessz(3,ib) = hessz(3,ib) + deddt*dzibzic
     &                               + d2eddt2*ddtdzic*ddtdzib
                  hessx(1,ia) = hessx(1,ia) + deddt*dxiaxic
     &                               + d2eddt2*ddtdxic*ddtdxia
                  hessx(2,ia) = hessx(2,ia) + deddt*dyiaxic
     &                               + d2eddt2*ddtdxic*ddtdyia
                  hessx(3,ia) = hessx(3,ia) + deddt*dziaxic
     &                               + d2eddt2*ddtdxic*ddtdzia
                  hessy(1,ia) = hessy(1,ia) + deddt*dxiayic
     &                               + d2eddt2*ddtdyic*ddtdxia
                  hessy(2,ia) = hessy(2,ia) + deddt*dyiayic
     &                               + d2eddt2*ddtdyic*ddtdyia
                  hessy(3,ia) = hessy(3,ia) + deddt*dziayic
     &                               + d2eddt2*ddtdyic*ddtdzia
                  hessz(1,ia) = hessz(1,ia) + deddt*dxiazic
     &                               + d2eddt2*ddtdzic*ddtdxia
                  hessz(2,ia) = hessz(2,ia) + deddt*dyiazic
     &                               + d2eddt2*ddtdzic*ddtdyia
                  hessz(3,ia) = hessz(3,ia) + deddt*dziazic
     &                               + d2eddt2*ddtdzic*ddtdzia
               end if
c
c     construct a second orthogonal direction for linear angles
c
               if (linear) then
                  linear = .false.
                  xpo = xp
                  ypo = yp
                  zpo = zp
                  xp = ypo*zab - zpo*yab
                  yp = zpo*xab - xpo*zab
                  zp = xpo*yab - ypo*xab
                  rp = sqrt(xp*xp + yp*yp + zp*zp)
                  goto 10
               end if
            end if
         end if
      end do
c
c     compute the Hessian elements for torsion restraints
c
      do ktors = 1, ntfix
         ia = itfix(1,ktors)
         ib = itfix(2,ktors)
         ic = itfix(3,ktors)
         id = itfix(4,ktors)
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,0,0)
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
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     calculate the pseudoenergy master chain rule terms
c
               force = tfix(1,ktors)
               tf1 = tfix(2,ktors)
               tf2 = tfix(3,ktors)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               dt = angle - target
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               dedphi = 2.0d0 * radian * force * dt
               d2edphi2 = 2.0d0 * radian**2 * force
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
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
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
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
c     chain rule terms for first derivative components
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
     &                              + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                              + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                              + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxia*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                              + d2edphi2*dphidyia*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                              + d2edphi2*dphidzia*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                              + d2edphi2*dphidxia*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                              + d2edphi2*dphidyia*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                              + d2edphi2*dphidzia*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                              + d2edphi2*dphidxia*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                              + d2edphi2*dphidyia*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                              + d2edphi2*dphidzia*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxia*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                              + d2edphi2*dphidyia*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                              + d2edphi2*dphidzia*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                              + d2edphi2*dphidxia*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                              + d2edphi2*dphidyia*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                              + d2edphi2*dphidzia*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                              + d2edphi2*dphidxia*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                              + d2edphi2*dphidyia*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                              + d2edphi2*dphidzia*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                              + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                              + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                              + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                              + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                              + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                              + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                              + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                              + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                              + d2edphi2*dphidxib*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                              + d2edphi2*dphidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                              + d2edphi2*dphidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                              + d2edphi2*dphidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                              + d2edphi2*dphidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                              + d2edphi2*dphidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                              + d2edphi2*dphidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                              + d2edphi2*dphidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                              + d2edphi2*dphidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                              + d2edphi2*dphidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                              + d2edphi2*dphidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                              + d2edphi2*dphidxib*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                              + d2edphi2*dphidyib*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                              + d2edphi2*dphidzib*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                              + d2edphi2*dphidxib*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                              + d2edphi2*dphidyib*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                              + d2edphi2*dphidzib*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                              + d2edphi2*dphidxib*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                              + d2edphi2*dphidyib*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                              + d2edphi2*dphidzib*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                              + d2edphi2*dphidxib*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                              + d2edphi2*dphidyib*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                              + d2edphi2*dphidzib*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                              + d2edphi2*dphidxib*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                              + d2edphi2*dphidyib*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                              + d2edphi2*dphidzib*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                              + d2edphi2*dphidxib*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                              + d2edphi2*dphidyib*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                              + d2edphi2*dphidzib*dphidzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                              + d2edphi2*dphidxic*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                              + d2edphi2*dphidyic*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                              + d2edphi2*dphidzic*dphidzic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxic*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                              + d2edphi2*dphidyic*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                              + d2edphi2*dphidzic*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                              + d2edphi2*dphidxic*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                              + d2edphi2*dphidyic*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                              + d2edphi2*dphidzic*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                              + d2edphi2*dphidxic*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                              + d2edphi2*dphidyic*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                              + d2edphi2*dphidzic*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                              + d2edphi2*dphidxic*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                              + d2edphi2*dphidyic*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                              + d2edphi2*dphidzic*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                              + d2edphi2*dphidxic*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                              + d2edphi2*dphidyic*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                              + d2edphi2*dphidzic*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                              + d2edphi2*dphidxic*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                              + d2edphi2*dphidyic*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                              + d2edphi2*dphidzic*dphidzib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                              + d2edphi2*dphidxic*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                              + d2edphi2*dphidyic*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                              + d2edphi2*dphidzic*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                              + d2edphi2*dphidxic*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                              + d2edphi2*dphidyic*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                              + d2edphi2*dphidzic*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                              + d2edphi2*dphidxic*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                              + d2edphi2*dphidyic*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                              + d2edphi2*dphidzic*dphidzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                              + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                              + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                              + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                              + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                              + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                              + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                              + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                              + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                              + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                              + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                              + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                              + d2edphi2*dphidxid*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                              + d2edphi2*dphidyid*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                              + d2edphi2*dphidzid*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                              + d2edphi2*dphidxid*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                              + d2edphi2*dphidyid*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                              + d2edphi2*dphidzid*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                              + d2edphi2*dphidxid*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                              + d2edphi2*dphidyid*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                              + d2edphi2*dphidzid*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                              + d2edphi2*dphidxid*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                              + d2edphi2*dphidyid*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                              + d2edphi2*dphidzid*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                              + d2edphi2*dphidxid*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                              + d2edphi2*dphidyid*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                              + d2edphi2*dphidzid*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                              + d2edphi2*dphidxid*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                              + d2edphi2*dphidyid*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                              + d2edphi2*dphidzid*dphidzic
               end if
            end if
         end if
      end do
c
c     compute the Hessian elements for group distance restraints
c
      do kdist = 1, ngfix
         ia = igfix(1,kdist)
         ib = igfix(2,kdist)
         proceed = (grplist(i).eq.ia .or. grplist(i).eq.ib)
         if (proceed) then
            if (grplist(i) .eq. ib) then
               ib = ia
               ia = grplist(i)
            end if
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            do j = igrp(1,ia), igrp(2,ia)
              k = kgrp(j)
              weigh = mass(k)
              xcm = xcm + x(k)*weigh
              ycm = ycm + y(k)*weigh
              zcm = zcm + z(k)*weigh
            end do
            weigha = max(1.0d0,grpmass(ia))
            xr = xcm / weigha
            yr = ycm / weigha
            zr = zcm / weigha
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            do j = igrp(1,ib), igrp(2,ib)
              k = kgrp(j)
              weigh = mass(k)
              xcm = xcm + x(k)*weigh
              ycm = ycm + y(k)*weigh
              zcm = zcm + z(k)*weigh
            end do
            weighb = max(1.0d0,grpmass(ib))
            xr = xr - xcm/weighb
            yr = yr - ycm/weighb
            zr = zr - zcm/weighb
            intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &                  molcule(kgrp(igrp(1,ib))))
            if (use_bounds .and. intermol)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            force = gfix(1,kdist)
            gf1 = gfix(2,kdist)
            gf2 = gfix(3,kdist)
            target = r
            if (r .lt. gf1)  target = gf1
            if (r .gt. gf2)  target = gf2
            dt = r - target
            deddt = 2.0d0 * force
c
c     set the chain rule terms for the Hessian elements
c
            if (r .eq. 0.0d0) then
               r = 0.0001d0
               r2 = r * r
            end if
            de = deddt * dt/r
            term = (deddt-de) / r2
            termx = term * xr
            termy = term * yr
            termz = term * zr
            d2e(1,1) = termx*xr + de
            d2e(1,2) = termx*yr
            d2e(1,3) = termx*zr
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yr + de
            d2e(2,3) = termy*zr
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do k = igrp(1,ia), igrp(2,ia)
               m = kgrp(k)
               ratio = mass(i)*mass(m) / (weigha*weigha)
               do j = 1, 3
                  hessx(j,m) = hessx(j,m) + d2e(1,j)*ratio
                  hessy(j,m) = hessy(j,m) + d2e(2,j)*ratio
                  hessz(j,m) = hessz(j,m) + d2e(3,j)*ratio
               end do
            end do
            do k = igrp(1,ib), igrp(2,ib)
               m = kgrp(k)
               ratio = mass(i)*mass(m) / (weigha*weighb)
               do j = 1, 3
                  hessx(j,m) = hessx(j,m) - d2e(1,j)*ratio
                  hessy(j,m) = hessy(j,m) - d2e(2,j)*ratio
                  hessz(j,m) = hessz(j,m) - d2e(3,j)*ratio
               end do
            end do
         end if
      end do
c
c     compute the Hessian elements for chirality restraints
c
      do kchir = 1, nchir
         ia = ichir(1,kchir)
         ib = ichir(2,kchir)
         ic = ichir(3,kchir)
         id = ichir(4,kchir)
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed) then
            xad = x(ia) - x(id)
            yad = y(ia) - y(id)
            zad = z(ia) - z(id)
            xbd = x(ib) - x(id)
            ybd = y(ib) - y(id)
            zbd = z(ib) - z(id)
            xcd = x(ic) - x(id)
            ycd = y(ic) - y(id)
            zcd = z(ic) - z(id)
            c1 = ybd*zcd - zbd*ycd
            c2 = ycd*zad - zcd*yad
            c3 = yad*zbd - zad*ybd
            vol = xad*c1 + xbd*c2 + xcd*c3
            force = chir(1,kchir)
            cf1 = chir(2,kchir)
            cf2 = chir(3,kchir)
            target = vol
            if (vol .lt. min(cf1,cf2))  target = min(cf1,cf2)
            if (vol .gt. max(cf1,cf2))  target = max(cf1,cf2)
            dt = vol - target
            dt2 = dt * dt
            deddt = 2.0d0 * force * dt
            d2eddt2 = 2.0d0 * force
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
               d2eddt2 = d2eddt2 * fgrp
            end if
c
c     chain rule terms for first derivative components
c
            term = sqrt(d2eddt2)
            ddtdxia = term * (ybd*zcd - zbd*ycd)
            ddtdyia = term * (zbd*xcd - xbd*zcd)
            ddtdzia = term * (xbd*ycd - ybd*xcd)
            ddtdxib = term * (zad*ycd - yad*zcd)
            ddtdyib = term * (xad*zcd - zad*xcd)
            ddtdzib = term * (yad*xcd - xad*ycd)
            ddtdxic = term * (yad*zbd - zad*ybd)
            ddtdyic = term * (zad*xbd - xad*zbd)
            ddtdzic = term * (xad*ybd - yad*xbd)
            ddtdxid = -ddtdxia - ddtdxib - ddtdxic
            ddtdyid = -ddtdyia - ddtdyib - ddtdyic
            ddtdzid = -ddtdzia - ddtdzib - ddtdzic
c
c     chain rule terms for second derivative components (*deddt)
c
            dyiaxib = -deddt * zcd
            dziaxib = deddt * ycd
            dxiayib = deddt * zcd
            dziayib = -deddt * xcd
            dxiazib = -deddt * ycd
            dyiazib = deddt * xcd
            dyiaxic = deddt * zbd
            dziaxic = -deddt * ybd
            dxiayic = -deddt * zbd
            dziayic = deddt * xbd
            dxiazic = deddt * ybd
            dyiazic = -deddt * xbd
            dyibxic = -deddt * zad
            dzibxic = deddt * yad
            dxibyic = deddt * zad
            dzibyic = -deddt * xad
            dxibzic = -deddt * yad
            dyibzic = deddt * xad
            dyiaxid = -dyiaxib - dyiaxic
            dziaxid = -dziaxib - dziaxic
            dxiayid = -dxiayib - dxiayic
            dziayid = -dziayib - dziayic
            dxiazid = -dxiazib - dxiazic
            dyiazid = -dyiazib - dyiazic
            dyibxid = -dxiayib - dyibxic
            dzibxid = -dxiazib - dzibxic
            dxibyid = -dyiaxib - dxibyic
            dzibyid = -dyiazib - dzibyic
            dxibzid = -dziaxib - dxibzic
            dyibzid = -dziayib - dyibzic
            dyicxid = -dxiayic - dxibyic
            dzicxid = -dxiazic - dxibzic
            dxicyid = -dyiaxic - dyibxic
            dzicyid = -dyiazic - dyibzic
            dxiczid = -dziaxic - dzibxic
            dyiczid = -dziayic - dzibyic
c
c     increment diagonal and off-diagonal Hessian elements
c
            if (i .eq. ia) then
               hessx(1,ia) = hessx(1,ia) + ddtdxia*ddtdxia
               hessy(1,ia) = hessy(1,ia) + ddtdxia*ddtdyia
               hessz(1,ia) = hessz(1,ia) + ddtdxia*ddtdzia
               hessx(2,ia) = hessx(2,ia) + ddtdxia*ddtdyia
               hessy(2,ia) = hessy(2,ia) + ddtdyia*ddtdyia
               hessz(2,ia) = hessz(2,ia) + ddtdyia*ddtdzia
               hessx(3,ia) = hessx(3,ia) + ddtdxia*ddtdzia
               hessy(3,ia) = hessy(3,ia) + ddtdyia*ddtdzia
               hessz(3,ia) = hessz(3,ia) + ddtdzia*ddtdzia
               hessx(1,ib) = hessx(1,ib) + ddtdxia*ddtdxib
               hessy(1,ib) = hessy(1,ib) + ddtdyia*ddtdxib + dyiaxib
               hessz(1,ib) = hessz(1,ib) + ddtdzia*ddtdxib + dziaxib
               hessx(2,ib) = hessx(2,ib) + ddtdxia*ddtdyib + dxiayib
               hessy(2,ib) = hessy(2,ib) + ddtdyia*ddtdyib
               hessz(2,ib) = hessz(2,ib) + ddtdzia*ddtdyib + dziayib
               hessx(3,ib) = hessx(3,ib) + ddtdxia*ddtdzib + dxiazib
               hessy(3,ib) = hessy(3,ib) + ddtdyia*ddtdzib + dyiazib
               hessz(3,ib) = hessz(3,ib) + ddtdzia*ddtdzib
               hessx(1,ic) = hessx(1,ic) + ddtdxia*ddtdxic
               hessy(1,ic) = hessy(1,ic) + ddtdyia*ddtdxic + dyiaxic
               hessz(1,ic) = hessz(1,ic) + ddtdzia*ddtdxic + dziaxic
               hessx(2,ic) = hessx(2,ic) + ddtdxia*ddtdyic + dxiayic
               hessy(2,ic) = hessy(2,ic) + ddtdyia*ddtdyic
               hessz(2,ic) = hessz(2,ic) + ddtdzia*ddtdyic + dziayic
               hessx(3,ic) = hessx(3,ic) + ddtdxia*ddtdzic + dxiazic
               hessy(3,ic) = hessy(3,ic) + ddtdyia*ddtdzic + dyiazic
               hessz(3,ic) = hessz(3,ic) + ddtdzia*ddtdzic
               hessx(1,id) = hessx(1,id) + ddtdxia*ddtdxid
               hessy(1,id) = hessy(1,id) + ddtdyia*ddtdxid + dyiaxid
               hessz(1,id) = hessz(1,id) + ddtdzia*ddtdxid + dziaxid
               hessx(2,id) = hessx(2,id) + ddtdxia*ddtdyid + dxiayid
               hessy(2,id) = hessy(2,id) + ddtdyia*ddtdyid
               hessz(2,id) = hessz(2,id) + ddtdzia*ddtdyid + dziayid
               hessx(3,id) = hessx(3,id) + ddtdxia*ddtdzid + dxiazid
               hessy(3,id) = hessy(3,id) + ddtdyia*ddtdzid + dyiazid
               hessz(3,id) = hessz(3,id) + ddtdzia*ddtdzid
            else if (i .eq. ib) then
               hessx(1,ib) = hessx(1,ib) + ddtdxib*ddtdxib
               hessy(1,ib) = hessy(1,ib) + ddtdxib*ddtdyib
               hessz(1,ib) = hessz(1,ib) + ddtdxib*ddtdzib
               hessx(2,ib) = hessx(2,ib) + ddtdxib*ddtdyib
               hessy(2,ib) = hessy(2,ib) + ddtdyib*ddtdyib
               hessz(2,ib) = hessz(2,ib) + ddtdyib*ddtdzib
               hessx(3,ib) = hessx(3,ib) + ddtdxib*ddtdzib
               hessy(3,ib) = hessy(3,ib) + ddtdyib*ddtdzib
               hessz(3,ib) = hessz(3,ib) + ddtdzib*ddtdzib
               hessx(1,ia) = hessx(1,ia) + ddtdxib*ddtdxia
               hessy(1,ia) = hessy(1,ia) + ddtdyib*ddtdxia + dxiayib
               hessz(1,ia) = hessz(1,ia) + ddtdzib*ddtdxia + dxiazib
               hessx(2,ia) = hessx(2,ia) + ddtdxib*ddtdyia + dyiaxib
               hessy(2,ia) = hessy(2,ia) + ddtdyib*ddtdyia
               hessz(2,ia) = hessz(2,ia) + ddtdzib*ddtdyia + dyiazib
               hessx(3,ia) = hessx(3,ia) + ddtdxib*ddtdzia + dziaxib
               hessy(3,ia) = hessy(3,ia) + ddtdyib*ddtdzia + dziayib
               hessz(3,ia) = hessz(3,ia) + ddtdzib*ddtdzia
               hessx(1,ic) = hessx(1,ic) + ddtdxib*ddtdxic
               hessy(1,ic) = hessy(1,ic) + ddtdyib*ddtdxic + dyibxic
               hessz(1,ic) = hessz(1,ic) + ddtdzib*ddtdxic + dzibxic
               hessx(2,ic) = hessx(2,ic) + ddtdxib*ddtdyic + dxibyic
               hessy(2,ic) = hessy(2,ic) + ddtdyib*ddtdyic
               hessz(2,ic) = hessz(2,ic) + ddtdzib*ddtdyic + dzibyic
               hessx(3,ic) = hessx(3,ic) + ddtdxib*ddtdzic + dxibzic
               hessy(3,ic) = hessy(3,ic) + ddtdyib*ddtdzic + dyibzic
               hessz(3,ic) = hessz(3,ic) + ddtdzib*ddtdzic
               hessx(1,id) = hessx(1,id) + ddtdxib*ddtdxid
               hessy(1,id) = hessy(1,id) + ddtdyib*ddtdxid + dyibxid
               hessz(1,id) = hessz(1,id) + ddtdzib*ddtdxid + dzibxid
               hessx(2,id) = hessx(2,id) + ddtdxib*ddtdyid + dxibyid
               hessy(2,id) = hessy(2,id) + ddtdyib*ddtdyid
               hessz(2,id) = hessz(2,id) + ddtdzib*ddtdyid + dzibyid
               hessx(3,id) = hessx(3,id) + ddtdxib*ddtdzid + dxibzid
               hessy(3,id) = hessy(3,id) + ddtdyib*ddtdzid + dyibzid
               hessz(3,id) = hessz(3,id) + ddtdzib*ddtdzid
            else if (i .eq. ic) then
               hessx(1,ic) = hessx(1,ic) + ddtdxic*ddtdxic
               hessy(1,ic) = hessy(1,ic) + ddtdxic*ddtdyic
               hessz(1,ic) = hessz(1,ic) + ddtdxic*ddtdzic
               hessx(2,ic) = hessx(2,ic) + ddtdxic*ddtdyic
               hessy(2,ic) = hessy(2,ic) + ddtdyic*ddtdyic
               hessz(2,ic) = hessz(2,ic) + ddtdyic*ddtdzic
               hessx(3,ic) = hessx(3,ic) + ddtdxic*ddtdzic
               hessy(3,ic) = hessy(3,ic) + ddtdyic*ddtdzic
               hessz(3,ic) = hessz(3,ic) + ddtdzic*ddtdzic
               hessx(1,ia) = hessx(1,ia) + ddtdxic*ddtdxia
               hessy(1,ia) = hessy(1,ia) + ddtdyic*ddtdxia + dxiayic
               hessz(1,ia) = hessz(1,ia) + ddtdzic*ddtdxia + dxiazic
               hessx(2,ia) = hessx(2,ia) + ddtdxic*ddtdyia + dyiaxic
               hessy(2,ia) = hessy(2,ia) + ddtdyic*ddtdyia
               hessz(2,ia) = hessz(2,ia) + ddtdzic*ddtdyia + dyiazic
               hessx(3,ia) = hessx(3,ia) + ddtdxic*ddtdzia + dziaxic
               hessy(3,ia) = hessy(3,ia) + ddtdyic*ddtdzia + dziayic
               hessz(3,ia) = hessz(3,ia) + ddtdzic*ddtdzia
               hessx(1,ib) = hessx(1,ib) + ddtdxic*ddtdxib
               hessy(1,ib) = hessy(1,ib) + ddtdyic*ddtdxib + dxibyic
               hessz(1,ib) = hessz(1,ib) + ddtdzic*ddtdxib + dxibzic
               hessx(2,ib) = hessx(2,ib) + ddtdxic*ddtdyib + dyibxic
               hessy(2,ib) = hessy(2,ib) + ddtdyic*ddtdyib
               hessz(2,ib) = hessz(2,ib) + ddtdzic*ddtdyib + dyibzic
               hessx(3,ib) = hessx(3,ib) + ddtdxic*ddtdzib + dzibxic
               hessy(3,ib) = hessy(3,ib) + ddtdyic*ddtdzib + dzibyic
               hessz(3,ib) = hessz(3,ib) + ddtdzic*ddtdzib
               hessx(1,id) = hessx(1,id) + ddtdxic*ddtdxid
               hessy(1,id) = hessy(1,id) + ddtdyic*ddtdxid + dyicxid
               hessz(1,id) = hessz(1,id) + ddtdzic*ddtdxid + dzicxid
               hessx(2,id) = hessx(2,id) + ddtdxic*ddtdyid + dxicyid
               hessy(2,id) = hessy(2,id) + ddtdyic*ddtdyid
               hessz(2,id) = hessz(2,id) + ddtdzic*ddtdyid + dzicyid
               hessx(3,id) = hessx(3,id) + ddtdxic*ddtdzid + dxiczid
               hessy(3,id) = hessy(3,id) + ddtdyic*ddtdzid + dyiczid
               hessz(3,id) = hessz(3,id) + ddtdzic*ddtdzid
            else if (i .eq. id) then
               hessx(1,id) = hessx(1,id) + ddtdxid*ddtdxid
               hessy(1,id) = hessy(1,id) + ddtdxid*ddtdyid
               hessz(1,id) = hessz(1,id) + ddtdxid*ddtdzid
               hessx(2,id) = hessx(2,id) + ddtdxid*ddtdyid
               hessy(2,id) = hessy(2,id) + ddtdyid*ddtdyid
               hessz(2,id) = hessz(2,id) + ddtdyid*ddtdzid
               hessx(3,id) = hessx(3,id) + ddtdxid*ddtdzid
               hessy(3,id) = hessy(3,id) + ddtdyid*ddtdzid
               hessz(3,id) = hessz(3,id) + ddtdzid*ddtdzid
               hessx(1,ia) = hessx(1,ia) + ddtdxid*ddtdxia
               hessy(1,ia) = hessy(1,ia) + ddtdyid*ddtdxia + dxiayid
               hessz(1,ia) = hessz(1,ia) + ddtdzid*ddtdxia + dxiazid
               hessx(2,ia) = hessx(2,ia) + ddtdxid*ddtdyia + dyiaxid
               hessy(2,ia) = hessy(2,ia) + ddtdyid*ddtdyia
               hessz(2,ia) = hessz(2,ia) + ddtdzid*ddtdyia + dyiazid
               hessx(3,ia) = hessx(3,ia) + ddtdxid*ddtdzia + dziaxid
               hessy(3,ia) = hessy(3,ia) + ddtdyid*ddtdzia + dziayid
               hessz(3,ia) = hessz(3,ia) + ddtdzid*ddtdzia
               hessx(1,ib) = hessx(1,ib) + ddtdxid*ddtdxib
               hessy(1,ib) = hessy(1,ib) + ddtdyid*ddtdxib + dxibyid
               hessz(1,ib) = hessz(1,ib) + ddtdzid*ddtdxib + dxibzid
               hessx(2,ib) = hessx(2,ib) + ddtdxid*ddtdyib + dyibxid
               hessy(2,ib) = hessy(2,ib) + ddtdyid*ddtdyib
               hessz(2,ib) = hessz(2,ib) + ddtdzid*ddtdyib + dyibzid
               hessx(3,ib) = hessx(3,ib) + ddtdxid*ddtdzib + dzibxid
               hessy(3,ib) = hessy(3,ib) + ddtdyid*ddtdzib + dzibyid
               hessz(3,ib) = hessz(3,ib) + ddtdzid*ddtdzib
               hessx(1,ic) = hessx(1,ic) + ddtdxid*ddtdxic
               hessy(1,ic) = hessy(1,ic) + ddtdyid*ddtdxic + dxicyid
               hessz(1,ic) = hessz(1,ic) + ddtdzid*ddtdxic + dxiczid
               hessx(2,ic) = hessx(2,ic) + ddtdxid*ddtdyic + dyicxid
               hessy(2,ic) = hessy(2,ic) + ddtdyid*ddtdyic
               hessz(2,ic) = hessz(2,ic) + ddtdzid*ddtdyic + dyiczid
               hessx(3,ic) = hessx(3,ic) + ddtdxid*ddtdzic + dzicxid
               hessy(3,ic) = hessy(3,ic) + ddtdyid*ddtdzic + dzicyid
               hessz(3,ic) = hessz(3,ic) + ddtdzid*ddtdzic
            end if
         end if
      end do
c
c     compute Hessian elements for a Gaussian basin restraint
c
      if (use_basin) then
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = 1, n
            proceed = (k .ne. i)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               term = -width * r2
               expterm = 0.0d0
               if (term .gt. -50.0d0)
     &            expterm = depth * width * exp(term)
               dedr = -2.0d0 * expterm
               d2edr2 = (-4.0d0*term-2.0d0) * expterm
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedr = dedr * fgrp
                  d2edr2 = d2edr2 * fgrp
               end if
c
c     set the chain rule terms for the Hessian elements
c
               if (r2 .eq. 0.0d0) then
                  term = 0.0d0
               else
                  term = (d2edr2-dedr) / r2
               end if
               termx = term * xr
               termy = term * yr
               termz = term * zr
               d2e(1,1) = termx*xr + dedr
               d2e(1,2) = termx*yr
               d2e(1,3) = termx*zr
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yr + dedr
               d2e(2,3) = termy*zr
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zr + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + d2e(1,j)
                  hessy(j,i) = hessy(j,i) + d2e(2,j)
                  hessz(j,i) = hessz(j,i) + d2e(3,j)
                  hessx(j,k) = hessx(j,k) - d2e(1,j)
                  hessy(j,k) = hessy(j,k) - d2e(2,j)
                  hessz(j,k) = hessz(j,k) - d2e(3,j)
               end do
            end if
         end do
      end if
c
c     compute Hessian elements for a spherical droplet restraint
c
      if (use_wall) then
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,0,0,0,0,0)
         if (proceed) then
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri2 = xi**2 + yi**2 + zi**2
            ri = sqrt(ri2)
            r = rwall + buffer - ri
            r2 = r * r
            r6 = r2 * r2 * r2
            r12 = r6 * r6
            if (ri .eq. 0.0d0) then
               ri = 1.0d0
               ri2 = 1.0d0
            end if
            dedr = (12.0d0*a/r12 - 6.0d0*b/r6) / (r*ri)
            d2edr2 = (156.0d0*a/r12 - 42.0d0*b/r6) / (r2*ri2)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               dedr = dedr * fgrp
               d2edr2 = d2edr2 * fgrp
            end if
c
c     set the chain rule terms for the Hessian elements
c
            d2edr2 = d2edr2 - dedr/ri2
            termx = d2edr2 * xi
            termy = d2edr2 * yi
            termz = d2edr2 * zi
            d2e(1,1) = termx*xi + dedr
            d2e(1,2) = termx*yi
            d2e(1,3) = termx*zi
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yi + dedr
            d2e(2,3) = termy*zi
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zi + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,i) = hessx(j,i) + d2e(1,j)
               hessy(j,i) = hessy(j,i) + d2e(2,j)
               hessz(j,i) = hessz(j,i) + d2e(3,j)
            end do
         end if
      end if
      return
      end
