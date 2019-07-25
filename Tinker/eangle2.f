c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangle2  --  atom-by-atom angle bend Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangle2" calculates second derivatives of the angle bending
c     energy for a single atom using a mixture of analytical and
c     finite difference methods; projected in-plane angles at trigonal
c     centers, special linear or Fourier angle bending terms are
c     optionally used
c
c
      subroutine eangle2 (i)
      use sizes
      use angbnd
      use angpot
      use atoms
      use group
      use hessn
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      real*8 eps,fgrp
      real*8 old,term
      real*8, allocatable :: de(:,:)
      real*8, allocatable :: d0(:,:)
      logical proceed
      logical twosided
c
c
c     compute analytical angle bending Hessian elements
c
      call eangle2a (i)
c
c     set stepsize for derivatives and default group weight
c
      eps = 1.0d-5
      fgrp = 1.0d0
      twosided = .false.
      if (n .le. 50)  twosided = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (de(3,n))
      allocate (d0(3,n))
c
c     calculate numerical in-plane bend Hessian for current atom
c
      do k = 1, nangle
         proceed = .false.
         if (angtyp(k) .eq. 'IN-PLANE') then
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
            if (proceed .and. use_group)
     &         call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         end if
         if (proceed) then
            term = fgrp / eps
c
c     find first derivatives for the base structure
c
            if (.not. twosided) then
               call eangle2b (k,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            if (twosided) then
               x(i) = x(i) - 0.5d0*eps
               call eangle2b (k,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            x(i) = x(i) + eps
            call eangle2b (k,de)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(de(j,id)-d0(j,id))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            if (twosided) then
               y(i) = y(i) - 0.5d0*eps
               call eangle2b (k,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            y(i) = y(i) + eps
            call eangle2b (k,de)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(de(j,id)-d0(j,id))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            if (twosided) then
               z(i) = z(i) - 0.5d0*eps
               call eangle2b (k,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            z(i) = z(i) + eps
            call eangle2b (k,de)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(de(j,id)-d0(j,id))
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (de)
      deallocate (d0)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eangle2a  --  angle bending Hessian; analytical  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eangle2a" calculates bond angle bending potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine eangle2a (iatom)
      use sizes
      use angbnd
      use angpot
      use atoms
      use bound
      use group
      use hessn
      use math
      implicit none
      integer i,iatom
      integer ia,ib,ic
      real*8 ideal,force
      real*8 fold,factor,dot
      real*8 cosine,sine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rcb2
      real*8 xpo,ypo,zpo
      real*8 xp,yp,zp,rp,rp2
      real*8 xrab,yrab,zrab
      real*8 xrcb,yrcb,zrcb
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 deddt,d2eddt2
      real*8 terma,termc
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 dxiaxia,dxiayia,dxiazia
      real*8 dxibxib,dxibyib,dxibzib
      real*8 dxicxic,dxicyic,dxiczic
      real*8 dyiayia,dyiazia,dziazia
      real*8 dyibyib,dyibzib,dzibzib
      real*8 dyicyic,dyiczic,dziczic
      real*8 dxibxia,dxibyia,dxibzia
      real*8 dyibxia,dyibyia,dyibzia
      real*8 dzibxia,dzibyia,dzibzia
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      logical proceed,linear
c
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ideal = anat(i)
         force = ak(i)
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
c     compute the bond angle bending Hessian elements
c
            if (angtyp(i) .ne. 'IN-PLANE') then
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
               if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
                  xp = ycb*zab - zcb*yab
                  yp = zcb*xab - xcb*zab
                  zp = xcb*yab - ycb*xab
                  rp = sqrt(xp*xp + yp*yp + zp*zp)
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     get the master chain rule terms for derivatives
c
                  if (angtyp(i) .eq. 'HARMONIC') then
                     dt = angle - ideal
                     dt2 = dt * dt
                     dt3 = dt2 * dt
                     dt4 = dt3 * dt
                     deddt = angunit * force * dt * radian
     &                         * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                             + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
                     d2eddt2 = angunit * force * radian**2
     &                        * (2.0d0 + 6.0d0*cang*dt + 12.0d0*qang*dt2
     &                            + 20.0d0*pang*dt3 + 30.0d0*sang*dt4)
                  else if (angtyp(i) .eq. 'LINEAR') then
                     factor = 2.0d0 * angunit * radian**2
                     sine = sqrt(1.0d0-cosine*cosine)
                     deddt = -factor * force * sine
                     d2eddt2 = -factor * force * cosine
                  else if (angtyp(i) .eq. 'FOURIER') then
                     fold = afld(i)
                     factor = 2.0d0 * angunit * (radian**2/fold)
                     cosine = cos((fold*angle-ideal)/radian)
                     sine = sin((fold*angle-ideal)/radian)
                     deddt = -factor * force * sine
                     d2eddt2 = -factor * force * fold * cosine
                  end if
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
                  if (rp .lt. 0.0001d0) then
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
   10             continue
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
                  if (ia .eq. iatom) then
                     hessx(1,ia) = hessx(1,ia) + deddt*dxiaxia
     &                                  + d2eddt2*ddtdxia*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dxiayia
     &                                  + d2eddt2*ddtdxia*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dxiazia
     &                                  + d2eddt2*ddtdxia*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dxiayia
     &                                  + d2eddt2*ddtdyia*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyiayia
     &                                  + d2eddt2*ddtdyia*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dyiazia
     &                                  + d2eddt2*ddtdyia*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dxiazia
     &                                  + d2eddt2*ddtdzia*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dyiazia
     &                                  + d2eddt2*ddtdzia*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dziazia
     &                                  + d2eddt2*ddtdzia*ddtdzia
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxia
     &                                  + d2eddt2*ddtdxia*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dyibxia
     &                                  + d2eddt2*ddtdxia*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dzibxia
     &                                  + d2eddt2*ddtdxia*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyia
     &                                  + d2eddt2*ddtdyia*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyia
     &                                  + d2eddt2*ddtdyia*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dzibyia
     &                                  + d2eddt2*ddtdyia*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzia
     &                                  + d2eddt2*ddtdzia*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzia
     &                                  + d2eddt2*ddtdzia*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzia
     &                                  + d2eddt2*ddtdzia*ddtdzib
                     hessx(1,ic) = hessx(1,ic) + deddt*dxiaxic
     &                                  + d2eddt2*ddtdxia*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxiayic
     &                                  + d2eddt2*ddtdxia*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxiazic
     &                                  + d2eddt2*ddtdxia*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dyiaxic
     &                                  + d2eddt2*ddtdyia*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyiayic
     &                                  + d2eddt2*ddtdyia*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyiazic
     &                                  + d2eddt2*ddtdyia*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dziaxic
     &                                  + d2eddt2*ddtdzia*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dziayic
     &                                  + d2eddt2*ddtdzia*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dziazic
     &                                  + d2eddt2*ddtdzia*ddtdzic
                  else if (ib .eq. iatom) then
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxib
     &                                  + d2eddt2*ddtdxib*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dxibyib
     &                                  + d2eddt2*ddtdxib*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dxibzib
     &                                  + d2eddt2*ddtdxib*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyib
     &                                  + d2eddt2*ddtdyib*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyib
     &                                  + d2eddt2*ddtdyib*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dyibzib
     &                                  + d2eddt2*ddtdyib*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzib
     &                                  + d2eddt2*ddtdzib*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzib
     &                                  + d2eddt2*ddtdzib*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzib
     &                                  + d2eddt2*ddtdzib*ddtdzib
                     hessx(1,ia) = hessx(1,ia) + deddt*dxibxia
     &                                  + d2eddt2*ddtdxib*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dxibyia
     &                                  + d2eddt2*ddtdxib*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dxibzia
     &                                  + d2eddt2*ddtdxib*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dyibxia
     &                                  + d2eddt2*ddtdyib*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyibyia
     &                                  + d2eddt2*ddtdyib*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dyibzia
     &                                  + d2eddt2*ddtdyib*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dzibxia
     &                                  + d2eddt2*ddtdzib*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dzibyia
     &                                  + d2eddt2*ddtdzib*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dzibzia
     &                                  + d2eddt2*ddtdzib*ddtdzia
                     hessx(1,ic) = hessx(1,ic) + deddt*dxibxic
     &                                  + d2eddt2*ddtdxib*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxibyic
     &                                  + d2eddt2*ddtdxib*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxibzic
     &                                  + d2eddt2*ddtdxib*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dyibxic
     &                                  + d2eddt2*ddtdyib*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyibyic
     &                                  + d2eddt2*ddtdyib*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyibzic
     &                                  + d2eddt2*ddtdyib*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dzibxic
     &                                  + d2eddt2*ddtdzib*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dzibyic
     &                                  + d2eddt2*ddtdzib*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dzibzic
     &                                  + d2eddt2*ddtdzib*ddtdzic
                  else if (ic .eq. iatom) then
                     hessx(1,ic) = hessx(1,ic) + deddt*dxicxic
     &                                  + d2eddt2*ddtdxic*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxicyic
     &                                  + d2eddt2*ddtdxic*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxiczic
     &                                  + d2eddt2*ddtdxic*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dxicyic
     &                                  + d2eddt2*ddtdyic*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyicyic
     &                                  + d2eddt2*ddtdyic*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyiczic
     &                                  + d2eddt2*ddtdyic*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dxiczic
     &                                  + d2eddt2*ddtdzic*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dyiczic
     &                                  + d2eddt2*ddtdzic*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dziczic
     &                                  + d2eddt2*ddtdzic*ddtdzic
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxic
     &                                  + d2eddt2*ddtdxic*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dyibxic
     &                                  + d2eddt2*ddtdxic*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dzibxic
     &                                  + d2eddt2*ddtdxic*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyic
     &                                  + d2eddt2*ddtdyic*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyic
     &                                  + d2eddt2*ddtdyic*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dzibyic
     &                                  + d2eddt2*ddtdyic*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzic
     &                                  + d2eddt2*ddtdzic*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzic
     &                                  + d2eddt2*ddtdzic*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzic
     &                                  + d2eddt2*ddtdzic*ddtdzib
                     hessx(1,ia) = hessx(1,ia) + deddt*dxiaxic
     &                                  + d2eddt2*ddtdxic*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dyiaxic
     &                                  + d2eddt2*ddtdxic*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dziaxic
     &                                  + d2eddt2*ddtdxic*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dxiayic
     &                                  + d2eddt2*ddtdyic*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyiayic
     &                                  + d2eddt2*ddtdyic*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dziayic
     &                                  + d2eddt2*ddtdyic*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dxiazic
     &                                  + d2eddt2*ddtdzic*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dyiazic
     &                                  + d2eddt2*ddtdzic*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dziazic
     &                                  + d2eddt2*ddtdzic*ddtdzia
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
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle2b  --  in-plane bend Hessian; numerical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle2b" computes projected in-plane bending first derivatives
c     for a single angle with respect to Cartesian coordinates;
c     used in computation of finite difference second derivatives
c
c
      subroutine eangle2b (i,de)
      use sizes
      use angbnd
      use angpot
      use atoms
      use bound
      use math
      implicit none
      integer i,ia,ib,ic,id
      real*8 ideal,force
      real*8 dot,cosine,angle
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,term
      real*8 terma,termc
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap
      real*8 xcp,ycp,zcp
      real*8 rap2,rcp2
      real*8 xt,yt,zt
      real*8 rt2,ptrt2
      real*8 xm,ym,zm,rm
      real*8 delta,delta2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxip,dedyip,dedzip
      real*8 dpdxia,dpdyia,dpdzia
      real*8 dpdxic,dpdyic,dpdzic
      real*8 de(3,*)
c
c
c     set the atom numbers and parameters for this angle
c
      ia = iang(1,i)
      ib = iang(2,i)
      ic = iang(3,i)
      id = iang(4,i)
      ideal = anat(i)
      force = ak(i)
c
c     get the coordinates of the atoms in the angle
c
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
c
c     zero out the first derivative components
c
      de(1,ia) = 0.0d0
      de(2,ia) = 0.0d0
      de(3,ia) = 0.0d0
      de(1,ib) = 0.0d0
      de(2,ib) = 0.0d0
      de(3,ib) = 0.0d0
      de(1,ic) = 0.0d0
      de(2,ic) = 0.0d0
      de(3,ic) = 0.0d0
      de(1,id) = 0.0d0
      de(2,id) = 0.0d0
      de(3,id) = 0.0d0
c
c     compute the projected in-plane angle gradient
c
      xad = xia - xid
      yad = yia - yid
      zad = zia - zid
      xbd = xib - xid
      ybd = yib - yid
      zbd = zib - zid
      xcd = xic - xid
      ycd = yic - yid
      zcd = zic - zid
      if (use_polymer) then
         call image (xad,yad,zad)
         call image (xbd,ybd,zbd)
         call image (xcd,ycd,zcd)
      end if
      xt = yad*zcd - zad*ycd
      yt = zad*xcd - xad*zcd
      zt = xad*ycd - yad*xcd
      rt2 = xt*xt + yt*yt + zt*zt
      delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
      xip = xib + xt*delta
      yip = yib + yt*delta
      zip = zib + zt*delta
      xap = xia - xip
      yap = yia - yip
      zap = zia - zip
      xcp = xic - xip
      ycp = yic - yip
      zcp = zic - zip
      if (use_polymer) then
         call image (xap,yap,zap)
         call image (xcp,ycp,zcp)
      end if
      rap2 = xap*xap + yap*yap + zap*zap
      rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
      if (rap2.ne.0.0d0 .and. rcp2.ne.0.0d0) then
         xm = ycp*zap - zcp*yap
         ym = zcp*xap - xcp*zap
         zm = xcp*yap - ycp*xap
         rm = sqrt(xm*xm + ym*ym + zm*zm)
         rm = max(rm,0.0001d0)
         dot = xap*xcp + yap*ycp + zap*zcp
         cosine = dot / sqrt(rap2*rcp2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
c
c     get the master chain rule term for derivatives
c
         dt = angle - ideal
         dt2 = dt * dt
         dt3 = dt2 * dt
         dt4 = dt2 * dt2
         deddt = angunit * force * dt * radian
     &             * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                  + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     chain rule terms for first derivative components
c
         terma = -deddt / (rap2*rm)
         termc = deddt / (rcp2*rm)
         dedxia = terma * (yap*zm-zap*ym)
         dedyia = terma * (zap*xm-xap*zm)
         dedzia = terma * (xap*ym-yap*xm)
         dedxic = termc * (ycp*zm-zcp*ym)
         dedyic = termc * (zcp*xm-xcp*zm)
         dedzic = termc * (xcp*ym-ycp*xm)
         dedxip = -dedxia - dedxic
         dedyip = -dedyia - dedyic
         dedzip = -dedzia - dedzic
c
c     chain rule components for the projection of the central atom
c
         delta2 = 2.0d0 * delta
         ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
         term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
         dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
         term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
         dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
         term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
         dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
         term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
         dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
         term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
         dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
         term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
         dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     compute derivative components for this interaction
c
         dedxia = dedxia + dpdxia
         dedyia = dedyia + dpdyia
         dedzia = dedzia + dpdzia
         dedxib = dedxip
         dedyib = dedyip
         dedzib = dedzip
         dedxic = dedxic + dpdxic
         dedyic = dedyic + dpdyic
         dedzic = dedzic + dpdzic
         dedxid = -dedxia - dedxib - dedxic
         dedyid = -dedyia - dedyib - dedyic
         dedzid = -dedzia - dedzib - dedzic
c
c     set the in-plane angle bending first derivatives
c
         de(1,ia) = dedxia
         de(2,ia) = dedyia
         de(3,ia) = dedzia
         de(1,ib) = dedxib
         de(2,ib) = dedyib
         de(3,ib) = dedzib
         de(1,ic) = dedxic
         de(2,ic) = dedyic
         de(3,ic) = dedzic
         de(1,id) = dedxid
         de(2,id) = dedyid
         de(3,id) = dedzid
      end if
      return
      end
