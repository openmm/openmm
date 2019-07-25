c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epitors2  --  pi-system torsion Hessian; numer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epitors2" calculates the second derivatives of the pi-system
c     torsion energy for a single atom using finite difference methods
c
c
      subroutine epitors2 (i)
      use sizes
      use angbnd
      use atoms
      use group
      use hessn
      use pitors
      use usage
      implicit none
      integer i,j,ipitors
      integer ia,ib,ic
      integer id,ie,ig
      real*8 eps,fgrp
      real*8 old,term
      real*8, allocatable :: de(:,:)
      real*8, allocatable :: d0(:,:)
      logical proceed
      logical twosided
c
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
c     calculate numerical pi-system torsion Hessian for current atom
c
      do ipitors = 1, npitors
         ia = ipit(1,ipitors)
         ib = ipit(2,ipitors)
         ic = ipit(3,ipitors)
         id = ipit(4,ipitors)
         ie = ipit(5,ipitors)
         ig = ipit(6,ipitors)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,ig)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic) .or.
     &                              use(id) .or. use(ie) .or. use(ig))
         if (proceed) then
            term = fgrp / eps
c
c     find first derivatives for the base structure
c
            if (.not. twosided) then
               call epitors2a (ipitors,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
                  d0(j,ie) = de(j,ie)
                  d0(j,ig) = de(j,ig)
               end do
            end if
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            if (twosided) then
               x(i) = x(i) - 0.5d0*eps
               call epitors2a (ipitors,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
                  d0(j,ie) = de(j,ie)
                  d0(j,ig) = de(j,ig)
               end do
            end if
            x(i) = x(i) + eps
            call epitors2a (ipitors,de)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(de(j,id)-d0(j,id))
               hessx(j,ie) = hessx(j,ie) + term*(de(j,ie)-d0(j,ie))
               hessx(j,ig) = hessx(j,ig) + term*(de(j,ig)-d0(j,ig))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            if (twosided) then
               y(i) = y(i) - 0.5d0*eps
               call epitors2a (ipitors,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
                  d0(j,ie) = de(j,ie)
                  d0(j,ig) = de(j,ig)
               end do
            end if
            y(i) = y(i) + eps
            call epitors2a (ipitors,de)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(de(j,id)-d0(j,id))
               hessy(j,ie) = hessy(j,ie) + term*(de(j,ie)-d0(j,ie))
               hessy(j,ig) = hessy(j,ig) + term*(de(j,ig)-d0(j,ig))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            if (twosided) then
               z(i) = z(i) - 0.5d0*eps
               call epitors2a (ipitors,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
                  d0(j,ie) = de(j,ie)
                  d0(j,ig) = de(j,ig)
               end do
            end if
            z(i) = z(i) + eps
            call epitors2a (ipitors,de)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(de(j,id)-d0(j,id))
               hessz(j,ie) = hessz(j,ie) + term*(de(j,ie)-d0(j,ie))
               hessz(j,ig) = hessz(j,ig) + term*(de(j,ig)-d0(j,ig))
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epitors2a  --  pi-system torsion derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epitors2a" calculates the pi-system torsion first derivatives;
c     used in computation of finite difference second derivatives
c
c
      subroutine epitors2a (i,de)
      use sizes
      use atoms
      use bound
      use deriv
      use pitors
      use torpot
      implicit none
      integer i,ia,ib,ic
      integer id,ie,ig
      real*8 dedphi,dphi2
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu
      real*8 rdc,rtru
      real*8 v2,c2,s2
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xig,yig,zig
      real*8 xip,yip,zip
      real*8 xiq,yiq,ziq
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xec,yec,zec
      real*8 xgc,ygc,zgc
      real*8 xcp,ycp,zcp
      real*8 xdc,ydc,zdc
      real*8 xqd,yqd,zqd
      real*8 xdp,ydp,zdp
      real*8 xqc,yqc,zqc
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxie,dedyie,dedzie
      real*8 dedxig,dedyig,dedzig
      real*8 dedxip,dedyip,dedzip
      real*8 dedxiq,dedyiq,dedziq
      real*8 de(3,*)
c
c
c     set the atom numbers for this pi-system torsion
c
      ia = ipit(1,i)
      ib = ipit(2,i)
      ic = ipit(3,i)
      id = ipit(4,i)
      ie = ipit(5,i)
      ig = ipit(6,i)
c
c     compute the value of the pi-system torsion angle
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
      xie = x(ie)
      yie = y(ie)
      zie = z(ie)
      xig = x(ig)
      yig = y(ig)
      zig = z(ig)
      xad = xia - xid
      yad = yia - yid
      zad = zia - zid
      xbd = xib - xid
      ybd = yib - yid
      zbd = zib - zid
      xec = xie - xic
      yec = yie - yic
      zec = zie - zic
      xgc = xig - xic
      ygc = yig - yic
      zgc = zig - zic
      if (use_polymer) then
         call image (xad,yad,zad)
         call image (xbd,ybd,zbd)
         call image (xec,yec,zec)
         call image (xgc,ygc,zgc)
      end if
      xip = yad*zbd - ybd*zad + xic
      yip = zad*xbd - zbd*xad + yic
      zip = xad*ybd - xbd*yad + zic
      xiq = yec*zgc - ygc*zec + xid
      yiq = zec*xgc - zgc*xec + yid
      ziq = xec*ygc - xgc*yec + zid
      xcp = xic - xip
      ycp = yic - yip
      zcp = zic - zip
      xdc = xid - xic
      ydc = yid - yic
      zdc = zid - zic
      xqd = xiq - xid
      yqd = yiq - yid
      zqd = ziq - zid
      if (use_polymer) then
         call image (xcp,ycp,zcp)
         call image (xdc,ydc,zdc)
         call image (xqd,yqd,zqd)
      end if
      xt = ycp*zdc - ydc*zcp
      yt = zcp*xdc - zdc*xcp
      zt = xcp*ydc - xdc*ycp
      xu = ydc*zqd - yqd*zdc
      yu = zdc*xqd - zqd*xdc
      zu = xdc*yqd - xqd*ydc
      xtu = yt*zu - yu*zt
      ytu = zt*xu - zu*xt
      ztu = xt*yu - xu*yt
      rt2 = xt*xt + yt*yt + zt*zt
      ru2 = xu*xu + yu*yu + zu*zu
      rtru = sqrt(rt2 * ru2)
      if (rtru .ne. 0.0d0) then
         rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
         cosine = (xt*xu + yt*yu + zt*zu) / rtru
         sine = (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru)
c
c     set the pi-system torsion parameters for this angle
c
         v2 = kpit(i)
         c2 = -1.0d0
         s2 = 0.0d0
c
c     compute the multiple angle trigonometry and the phase terms
c
         cosine2 = cosine*cosine - sine*sine
         sine2 = 2.0d0 * cosine * sine
         dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
c
c     calculate pi-system torsion energy and master chain rule term
c
         dedphi = ptorunit * v2 * dphi2
c
c     chain rule terms for first derivative components
c
         xdp = xid - xip
         ydp = yid - yip
         zdp = zid - zip
         xqc = xiq - xic
         yqc = yiq - yic
         zqc = ziq - zic
         dedxt = dedphi * (yt*zdc - ydc*zt) / (rt2*rdc)
         dedyt = dedphi * (zt*xdc - zdc*xt) / (rt2*rdc)
         dedzt = dedphi * (xt*ydc - xdc*yt) / (rt2*rdc)
         dedxu = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc)
         dedyu = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc)
         dedzu = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc)
c
c     compute first derivative components for pi-system angle
c
         dedxip = zdc*dedyt - ydc*dedzt
         dedyip = xdc*dedzt - zdc*dedxt
         dedzip = ydc*dedxt - xdc*dedyt
         dedxic = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu
         dedyic = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu
         dedzic = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu
         dedxid = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu
         dedyid = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu
         dedzid = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu
         dedxiq = zdc*dedyu - ydc*dedzu
         dedyiq = xdc*dedzu - zdc*dedxu
         dedziq = ydc*dedxu - xdc*dedyu
c
c     compute first derivative components for individual atoms
c
         dedxia = ybd*dedzip - zbd*dedyip
         dedyia = zbd*dedxip - xbd*dedzip
         dedzia = xbd*dedyip - ybd*dedxip
         dedxib = zad*dedyip - yad*dedzip
         dedyib = xad*dedzip - zad*dedxip
         dedzib = yad*dedxip - xad*dedyip
         dedxie = ygc*dedziq - zgc*dedyiq
         dedyie = zgc*dedxiq - xgc*dedziq
         dedzie = xgc*dedyiq - ygc*dedxiq
         dedxig = zec*dedyiq - yec*dedziq
         dedyig = xec*dedziq - zec*dedxiq
         dedzig = yec*dedxiq - xec*dedyiq
         dedxic = dedxic + dedxip - dedxie - dedxig
         dedyic = dedyic + dedyip - dedyie - dedyig
         dedzic = dedzic + dedzip - dedzie - dedzig
         dedxid = dedxid + dedxiq - dedxia - dedxib
         dedyid = dedyid + dedyiq - dedyia - dedyib
         dedzid = dedzid + dedziq - dedzia - dedzib
c
c     increment the total pi-system torsion energy and gradient
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
         de(1,ie) = dedxie
         de(2,ie) = dedyie
         de(3,ie) = dedzie
         de(1,ig) = dedxig
         de(2,ig) = dedyig
         de(3,ig) = dedzig
      end if
      return
      end
