c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrbnd2  --  stretch-bend Hessian; analytical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrbnd2" calculates the stretch-bend potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd2 (iatom)
      use sizes
      use angbnd
      use angpot
      use atoms
      use bndstr
      use bound
      use group
      use hessn
      use math
      use strbnd
      implicit none
      integer i,j,k,iatom
      integer ia,ib,ic,istrbnd
      real*8 angle,fgrp
      real*8 dot,cosine
      real*8 force1,force2
      real*8 dt,dr,dr1,dr2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab,rcb,rab2,rcb2
      real*8 xp,yp,zp,rp,rp2
      real*8 term,term1,term2
      real*8 xrab,yrab,zrab
      real*8 xrcb,yrcb,zrcb
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxib,ddrdyib,ddrdzib
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dtxiaxia,dtxiayia,dtxiazia
      real*8 dtxibxib,dtxibyib,dtxibzib
      real*8 dtxicxic,dtxicyic,dtxiczic
      real*8 dtyiayia,dtyiazia,dtziazia
      real*8 dtyibyib,dtyibzib,dtzibzib
      real*8 dtyicyic,dtyiczic,dtziczic
      real*8 dtxibxia,dtxibyia,dtxibzia
      real*8 dtyibxia,dtyibyia,dtyibzia
      real*8 dtzibxia,dtzibyia,dtzibzia
      real*8 dtxibxic,dtxibyic,dtxibzic
      real*8 dtyibxic,dtyibyic,dtyibzic
      real*8 dtzibxic,dtzibyic,dtzibzic
      real*8 dtxiaxic,dtxiayic,dtxiazic
      real*8 dtyiaxic,dtyiayic,dtyiazic
      real*8 dtziaxic,dtziayic,dtziazic
      real*8 drxiaxia,drxiayia,drxiazia
      real*8 drxibxib,drxibyib,drxibzib
      real*8 drxicxic,drxicyic,drxiczic
      real*8 dryiayia,dryiazia,drziazia
      real*8 dryibyib,dryibzib,drzibzib
      real*8 dryicyic,dryiczic,drziczic
      logical proceed
c
c
c     calculate the stretch-bend interaction Hessian elements
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = (iatom.eq.ia .or. iatom.eq.ib .or. iatom.eq.ic)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,0,0,0)
c
c     get the coordinates of the atoms in the angle
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (min(rab2,rcb2) .ne. 0.0d0) then
               rab = sqrt(rab2)
               rcb = sqrt(rcb2)
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.0001d0)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     first derivatives of angle with respect to coordinates
c
               dt = angle - anat(i)
               term1 = -radian / (rab2*rp)
               term2 = radian / (rcb2*rp)
               ddtdxia = term1 * (yab*zp-zab*yp)
               ddtdyia = term1 * (zab*xp-xab*zp)
               ddtdzia = term1 * (xab*yp-yab*xp)
               ddtdxic = term2 * (ycb*zp-zcb*yp)
               ddtdyic = term2 * (zcb*xp-xcb*zp)
               ddtdzic = term2 * (xcb*yp-ycb*xp)
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
c     second derivatives of angle with respect to coordinates
c
               dtxiaxia = term1*(xab*xcb-dot) + ddtdxia*(xcbp-xrab)
               dtxiayia = term1*(zp+yab*xcb) + ddtdxia*(ycbp-yrab)
               dtxiazia = term1*(zab*xcb-yp) + ddtdxia*(zcbp-zrab)
               dtyiayia = term1*(yab*ycb-dot) + ddtdyia*(ycbp-yrab)
               dtyiazia = term1*(xp+zab*ycb) + ddtdyia*(zcbp-zrab)
               dtziazia = term1*(zab*zcb-dot) + ddtdzia*(zcbp-zrab)
               dtxicxic = term2*(dot-xab*xcb) - ddtdxic*(xabp+xrcb)
               dtxicyic = term2*(zp-ycb*xab) - ddtdxic*(yabp+yrcb)
               dtxiczic = -term2*(yp+zcb*xab) - ddtdxic*(zabp+zrcb)
               dtyicyic = term2*(dot-yab*ycb) - ddtdyic*(yabp+yrcb)
               dtyiczic = term2*(xp-zcb*yab) - ddtdyic*(zabp+zrcb)
               dtziczic = term2*(dot-zab*zcb) - ddtdzic*(zabp+zrcb)
               dtxiaxic = term1*(yab*yab+zab*zab) - ddtdxia*xabp
               dtxiayic = -term1*xab*yab - ddtdxia*yabp
               dtxiazic = -term1*xab*zab - ddtdxia*zabp
               dtyiaxic = -term1*xab*yab - ddtdyia*xabp
               dtyiayic = term1*(xab*xab+zab*zab) - ddtdyia*yabp
               dtyiazic = -term1*yab*zab - ddtdyia*zabp
               dtziaxic = -term1*xab*zab - ddtdzia*xabp
               dtziayic = -term1*yab*zab - ddtdzia*yabp
               dtziazic = term1*(xab*xab+yab*yab) - ddtdzia*zabp
c
c     more angle deviation derivatives resulting from symmetry
c
               dtxibxia = -dtxiaxia - dtxiaxic
               dtxibyia = -dtxiayia - dtyiaxic
               dtxibzia = -dtxiazia - dtziaxic
               dtyibxia = -dtxiayia - dtxiayic
               dtyibyia = -dtyiayia - dtyiayic
               dtyibzia = -dtyiazia - dtziayic
               dtzibxia = -dtxiazia - dtxiazic
               dtzibyia = -dtyiazia - dtyiazic
               dtzibzia = -dtziazia - dtziazic
               dtxibxic = -dtxicxic - dtxiaxic
               dtxibyic = -dtxicyic - dtxiayic
               dtxibzic = -dtxiczic - dtxiazic
               dtyibxic = -dtxicyic - dtyiaxic
               dtyibyic = -dtyicyic - dtyiayic
               dtyibzic = -dtyiczic - dtyiazic
               dtzibxic = -dtxiczic - dtziaxic
               dtzibyic = -dtyiczic - dtziayic
               dtzibzic = -dtziczic - dtziazic
               dtxibxib = -dtxibxia - dtxibxic
               dtxibyib = -dtxibyia - dtxibyic
               dtxibzib = -dtxibzia - dtxibzic
               dtyibyib = -dtyibyia - dtyibyic
               dtyibzib = -dtyibzia - dtyibzic
               dtzibzib = -dtzibzia - dtzibzic
c
c     compute the values of the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               term = stbnunit * force1
               dr1 = term*(rab-bl(j))
               term1 = term / rab
               term = stbnunit * force2
               dr2 = term*(rcb-bl(k))
               term2 = term / rcb
               dr = dr1 + dr2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dr = dr * fgrp
                  dr1 = dr1 * fgrp
                  dr2 = dr2 * fgrp
                  term1 = term1 * fgrp
                  term2 = term2 * fgrp
               end if
c
c     first derivatives of bond length with respect to coordinates
c
               ddrdxia = term1 * xab
               ddrdyia = term1 * yab
               ddrdzia = term1 * zab
               ddrdxic = term2 * xcb
               ddrdyic = term2 * ycb
               ddrdzic = term2 * zcb
               ddrdxib = -ddrdxia - ddrdxic
               ddrdyib = -ddrdyia - ddrdyic
               ddrdzib = -ddrdzia - ddrdzic
c
c     abbreviations used in defining chain rule terms
c
               xab = xab / rab
               yab = yab / rab
               zab = zab / rab
               xcb = xcb / rcb
               ycb = ycb / rcb
               zcb = zcb / rcb
c
c     second derivatives of bond length with respect to coordinates
c
               drxiaxia = term1 * (1.0d0-xab*xab)
               drxiayia = -term1 * xab*yab
               drxiazia = -term1 * xab*zab
               dryiayia = term1 * (1.0d0-yab*yab)
               dryiazia = -term1 * yab*zab
               drziazia = term1 * (1.0d0-zab*zab)
               drxicxic = term2 * (1.0d0-xcb*xcb)
               drxicyic = -term2 * xcb*ycb
               drxiczic = -term2 * xcb*zcb
               dryicyic = term2 * (1.0d0-ycb*ycb)
               dryiczic = -term2 * ycb*zcb
               drziczic = term2 * (1.0d0-zcb*zcb)
               drxibxib = drxiaxia + drxicxic
               drxibyib = drxiayia + drxicyic
               drxibzib = drxiazia + drxiczic
               dryibyib = dryiayia + dryicyic
               dryibzib = dryiazia + dryiczic
               drzibzib = drziazia + drziczic
c
c     increment diagonal and non-diagonal Hessian elements
c
               if (ia .eq. iatom) then
                  hessx(1,ia) = hessx(1,ia) + dt*drxiaxia + dr*dtxiaxia
     &                             + 2.0d0*ddtdxia*ddrdxia
                  hessx(2,ia) = hessx(2,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdxia*ddrdyia + ddtdyia*ddrdxia
                  hessx(3,ia) = hessx(3,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdxia*ddrdzia + ddtdzia*ddrdxia
                  hessy(1,ia) = hessy(1,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdyia*ddrdxia + ddtdxia*ddrdyia
                  hessy(2,ia) = hessy(2,ia) + dt*dryiayia + dr*dtyiayia
     &                             + 2.0d0*ddtdyia*ddrdyia
                  hessy(3,ia) = hessy(3,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdyia*ddrdzia + ddtdzia*ddrdyia
                  hessz(1,ia) = hessz(1,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdzia*ddrdxia + ddtdxia*ddrdzia
                  hessz(2,ia) = hessz(2,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdzia*ddrdyia + ddtdyia*ddrdzia
                  hessz(3,ia) = hessz(3,ia) + dt*drziazia + dr*dtziazia
     &                             + 2.0d0*ddtdzia*ddrdzia
                  hessx(1,ib) = hessx(1,ib) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxia*ddrdxib + ddtdxib*ddrdxia
                  hessx(2,ib) = hessx(2,ib) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxia*ddrdyib + ddtdyib*ddrdxia
                  hessx(3,ib) = hessx(3,ib) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxia*ddrdzib + ddtdzib*ddrdxia
                  hessy(1,ib) = hessy(1,ib) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyia*ddrdxib + ddtdxib*ddrdyia
                  hessy(2,ib) = hessy(2,ib) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyia*ddrdyib + ddtdyib*ddrdyia
                  hessy(3,ib) = hessy(3,ib) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyia*ddrdzib + ddtdzib*ddrdyia
                  hessz(1,ib) = hessz(1,ib) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzia*ddrdxib + ddtdxib*ddrdzia
                  hessz(2,ib) = hessz(2,ib) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzia*ddrdyib + ddtdyib*ddrdzia
                  hessz(3,ib) = hessz(3,ib) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzia*ddrdzib + ddtdzib*ddrdzia
                  hessx(1,ic) = hessx(1,ic) + dr*dtxiaxic
     &                             + ddtdxia*ddrdxic + ddtdxic*ddrdxia
                  hessx(2,ic) = hessx(2,ic) + dr*dtxiayic
     &                             + ddtdxia*ddrdyic + ddtdyic*ddrdxia
                  hessx(3,ic) = hessx(3,ic) + dr*dtxiazic
     &                             + ddtdxia*ddrdzic + ddtdzic*ddrdxia
                  hessy(1,ic) = hessy(1,ic) + dr*dtyiaxic
     &                             + ddtdyia*ddrdxic + ddtdxic*ddrdyia
                  hessy(2,ic) = hessy(2,ic) + dr*dtyiayic
     &                             + ddtdyia*ddrdyic + ddtdyic*ddrdyia
                  hessy(3,ic) = hessy(3,ic) + dr*dtyiazic
     &                             + ddtdyia*ddrdzic + ddtdzic*ddrdyia
                  hessz(1,ic) = hessz(1,ic) + dr*dtziaxic
     &                             + ddtdzia*ddrdxic + ddtdxic*ddrdzia
                  hessz(2,ic) = hessz(2,ic) + dr*dtziayic
     &                             + ddtdzia*ddrdyic + ddtdyic*ddrdzia
                  hessz(3,ic) = hessz(3,ic) + dr*dtziazic
     &                             + ddtdzia*ddrdzic + ddtdzic*ddrdzia
               else if (ib .eq. iatom) then
                  hessx(1,ib) = hessx(1,ib) + dt*drxibxib + dr*dtxibxib
     &                             + 2.0d0*ddtdxib*ddrdxib
                  hessx(2,ib) = hessx(2,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdxib*ddrdyib + ddtdyib*ddrdxib
                  hessx(3,ib) = hessx(3,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdxib*ddrdzib + ddtdzib*ddrdxib
                  hessy(1,ib) = hessy(1,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdyib*ddrdxib + ddtdxib*ddrdyib
                  hessy(2,ib) = hessy(2,ib) + dt*dryibyib + dr*dtyibyib
     &                             + 2.0d0*ddtdyib*ddrdyib
                  hessy(3,ib) = hessy(3,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdyib*ddrdzib + ddtdzib*ddrdyib
                  hessz(1,ib) = hessz(1,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdzib*ddrdxib + ddtdxib*ddrdzib
                  hessz(2,ib) = hessz(2,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdzib*ddrdyib + ddtdyib*ddrdzib
                  hessz(3,ib) = hessz(3,ib) + dt*drzibzib + dr*dtzibzib
     &                             + 2.0d0*ddtdzib*ddrdzib
                  hessx(1,ia) = hessx(1,ia) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxib*ddrdxia + ddtdxia*ddrdxib
                  hessx(2,ia) = hessx(2,ia) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxib*ddrdyia + ddtdyia*ddrdxib
                  hessx(3,ia) = hessx(3,ia) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxib*ddrdzia + ddtdzia*ddrdxib
                  hessy(1,ia) = hessy(1,ia) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyib*ddrdxia + ddtdxia*ddrdyib
                  hessy(2,ia) = hessy(2,ia) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyib*ddrdyia + ddtdyia*ddrdyib
                  hessy(3,ia) = hessy(3,ia) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyib*ddrdzia + ddtdzia*ddrdyib
                  hessz(1,ia) = hessz(1,ia) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzib*ddrdxia + ddtdxia*ddrdzib
                  hessz(2,ia) = hessz(2,ia) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzib*ddrdyia + ddtdyia*ddrdzib
                  hessz(3,ia) = hessz(3,ia) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzib*ddrdzia + ddtdzia*ddrdzib
                  hessx(1,ic) = hessx(1,ic) - dt*drxicxic + dr*dtxibxic
     &                             + ddtdxib*ddrdxic + ddtdxic*ddrdxib
                  hessx(2,ic) = hessx(2,ic) - dt*drxicyic + dr*dtxibyic
     &                             + ddtdxib*ddrdyic + ddtdyic*ddrdxib
                  hessx(3,ic) = hessx(3,ic) - dt*drxiczic + dr*dtxibzic
     &                             + ddtdxib*ddrdzic + ddtdzic*ddrdxib
                  hessy(1,ic) = hessy(1,ic) - dt*drxicyic + dr*dtyibxic
     &                             + ddtdyib*ddrdxic + ddtdxic*ddrdyib
                  hessy(2,ic) = hessy(2,ic) - dt*dryicyic + dr*dtyibyic
     &                             + ddtdyib*ddrdyic + ddtdyic*ddrdyib
                  hessy(3,ic) = hessy(3,ic) - dt*dryiczic + dr*dtyibzic
     &                             + ddtdyib*ddrdzic + ddtdzic*ddrdyib
                  hessz(1,ic) = hessz(1,ic) - dt*drxiczic + dr*dtzibxic
     &                             + ddtdzib*ddrdxic + ddtdxic*ddrdzib
                  hessz(2,ic) = hessz(2,ic) - dt*dryiczic + dr*dtzibyic
     &                             + ddtdzib*ddrdyic + ddtdyic*ddrdzib
                  hessz(3,ic) = hessz(3,ic) - dt*drziczic + dr*dtzibzic
     &                             + ddtdzib*ddrdzic + ddtdzic*ddrdzib
               else if (ic .eq. iatom) then
                  hessx(1,ic) = hessx(1,ic) + dt*drxicxic + dr*dtxicxic
     &                             + 2.0d0*ddtdxic*ddrdxic
                  hessx(2,ic) = hessx(2,ic) + dt*drxicyic + dr*dtxicyic
     &                             + ddtdxic*ddrdyic + ddtdyic*ddrdxic
                  hessx(3,ic) = hessx(3,ic) + dt*drxiczic + dr*dtxiczic
     &                             + ddtdxic*ddrdzic + ddtdzic*ddrdxic
                  hessy(1,ic) = hessy(1,ic) + dt*drxicyic + dr*dtxicyic
     &                             + ddtdyic*ddrdxic + ddtdxic*ddrdyic
                  hessy(2,ic) = hessy(2,ic) + dt*dryicyic + dr*dtyicyic
     &                             + 2.0d0*ddtdyic*ddrdyic
                  hessy(3,ic) = hessy(3,ic) + dt*dryiczic + dr*dtyiczic
     &                             + ddtdyic*ddrdzic + ddtdzic*ddrdyic
                  hessz(1,ic) = hessz(1,ic) + dt*drxiczic + dr*dtxiczic
     &                             + ddtdzic*ddrdxic + ddtdxic*ddrdzic
                  hessz(2,ic) = hessz(2,ic) + dt*dryiczic + dr*dtyiczic
     &                             + ddtdzic*ddrdyic + ddtdyic*ddrdzic
                  hessz(3,ic) = hessz(3,ic) + dt*drziczic + dr*dtziczic
     &                             + 2.0d0*ddtdzic*ddrdzic
                  hessx(1,ib) = hessx(1,ib) - dt*drxicxic + dr*dtxibxic
     &                             + ddtdxic*ddrdxib + ddtdxib*ddrdxic
                  hessx(2,ib) = hessx(2,ib) - dt*drxicyic + dr*dtxibyic
     &                             + ddtdxic*ddrdyib + ddtdyib*ddrdxic
                  hessx(3,ib) = hessx(3,ib) - dt*drxiczic + dr*dtxibzic
     &                             + ddtdxic*ddrdzib + ddtdzib*ddrdxic
                  hessy(1,ib) = hessy(1,ib) - dt*drxicyic + dr*dtyibxic
     &                             + ddtdyic*ddrdxib + ddtdxib*ddrdyic
                  hessy(2,ib) = hessy(2,ib) - dt*dryicyic + dr*dtyibyic
     &                             + ddtdyic*ddrdyib + ddtdyib*ddrdyic
                  hessy(3,ib) = hessy(3,ib) - dt*dryiczic + dr*dtyibzic
     &                             + ddtdyic*ddrdzib + ddtdzib*ddrdyic
                  hessz(1,ib) = hessz(1,ib) - dt*drxiczic + dr*dtzibxic
     &                             + ddtdzic*ddrdxib + ddtdxib*ddrdzic
                  hessz(2,ib) = hessz(2,ib) - dt*dryiczic + dr*dtzibyic
     &                             + ddtdzic*ddrdyib + ddtdyib*ddrdzic
                  hessz(3,ib) = hessz(3,ib) - dt*drziczic + dr*dtzibzic
     &                             + ddtdzic*ddrdzib + ddtdzib*ddrdzic
                  hessx(1,ia) = hessx(1,ia) + dr*dtxiaxic
     &                             + ddtdxic*ddrdxia + ddtdxia*ddrdxic
                  hessx(2,ia) = hessx(2,ia) + dr*dtyiaxic
     &                             + ddtdxic*ddrdyia + ddtdyia*ddrdxic
                  hessx(3,ia) = hessx(3,ia) + dr*dtziaxic
     &                             + ddtdxic*ddrdzia + ddtdzia*ddrdxic
                  hessy(1,ia) = hessy(1,ia) + dr*dtxiayic
     &                             + ddtdyic*ddrdxia + ddtdxia*ddrdyic
                  hessy(2,ia) = hessy(2,ia) + dr*dtyiayic
     &                             + ddtdyic*ddrdyia + ddtdyia*ddrdyic
                  hessy(3,ia) = hessy(3,ia) + dr*dtziayic
     &                             + ddtdyic*ddrdzia + ddtdzia*ddrdyic
                  hessz(1,ia) = hessz(1,ia) + dr*dtxiazic
     &                             + ddtdzic*ddrdxia + ddtdxia*ddrdzic
                  hessz(2,ia) = hessz(2,ia) + dr*dtyiazic
     &                             + ddtdzic*ddrdyia + ddtdyia*ddrdzic
                  hessz(3,ia) = hessz(3,ia) + dr*dtziazic
     &                             + ddtdzic*ddrdzia + ddtdzia*ddrdzic
               end if
            end if
         end if
      end do
      return
      end
