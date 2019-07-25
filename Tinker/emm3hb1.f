c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb1  --  MM3 vdw & hbond energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb1" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding energy with respect to Cartesian
c     coordinates
c
c     literature references:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. II", Journal of Computational Chemistry,
c     19, 1001-1016 (1998)
c
c
      subroutine emm3hb1
      use energi
      use limits
      use vdwpot
      use virial
      implicit none
      real*8 elrc,vlrc
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_lights) then
         call emm3hb1b
      else if (use_vlist) then
         call emm3hb1c
      else
         call emm3hb1a
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr1 (elrc,vlrc)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb1a  --  double loop MM3 vdw-hbond derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb1a" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding energy with respect to Cartesian
c     coordinates using a pairwise double loop
c
c
      subroutine emm3hb1a
      use sizes
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer, allocatable :: iv14(:)
      real*8 e,de,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rvterm
      real*8 taper,dtaper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 dot,cosine,sine
      real*8 fterm,fcbuck,term
      real*8 deddr,ideal,ratio
      real*8 deddt,terma,termc
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei,use_hb
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     special cutoffs for very short and very long range terms
c
      expmin2 = 0.01d0
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (abuck*exp(-bbuck/expcut) - cbuck*(expcut**6))
     &                               / (expcut**12)
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
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
                     call image (xcb,ycb,zcb)
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                     deddt = -sine * ratio
                     deddr = cosine / (rab*ideal)
                  end if
                  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                     de = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bbuck / rv
                     expterm = abuck * exp(-bbuck/p)
                     fcbuck = fterm * cbuck * p6
                     e = eps * (expterm - fcbuck)
                     de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (use_hb) then
                        deddt = deddt * fgrp
                        deddr = deddr * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     find the chain rule terms for hydrogen bonding components
c
                  if (use_hb) then
                     term = eps * cbuck * p6
                     deddt = deddt * term
                     deddr = deddr * term
                     if (rik2 .gt. cut2) then
                        deddt = deddt * taper
                        deddr = deddr * taper
                     end if
                     terma = deddt / (rab2*max(rp,1.0d-6))
                     termc = -deddt / (rcb2*max(rp,1.0d-6))
                     dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                     dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                     dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                     dedxic = termc * (ycb*zp-zcb*yp)
                     dedyic = termc * (zcb*xp-xcb*zp)
                     dedzic = termc * (xcb*yp-ycb*xp)
                     dedxib = -dedxia - dedxic
                     dedyib = -dedyia - dedyic
                     dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                     dev(1,ia) = dev(1,ia) + dedxia
                     dev(2,ia) = dev(2,ia) + dedyia
                     dev(3,ia) = dev(3,ia) + dedzia
                     dev(1,ib) = dev(1,ib) + dedxib
                     dev(2,ib) = dev(2,ib) + dedyib
                     dev(3,ib) = dev(3,ib) + dedzib
                     dev(1,ic) = dev(1,ic) + dedxic
                     dev(2,ic) = dev(2,ic) + dedyic
                     dev(3,ic) = dev(3,ic) + dedzic
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  if (use_hb) then
                     vxx = xab*dedxia + xcb*dedxic
                     vyx = yab*dedxia + ycb*dedxic
                     vzx = zab*dedxia + zcb*dedxic
                     vyy = yab*dedyia + ycb*dedyic
                     vzy = zab*dedyia + zcb*dedyic
                     vzz = zab*dedzia + zcb*dedzic
                  end if
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               do j = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call imager (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     use_hb = .false.
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
                        use_hb = .true.
                        rv = radhbnd(kt,it)
                        eps = epshbnd(kt,it) / dielec
                        if (atomic(i) .eq. 1) then
                           ia = i
                           ib = i12(1,i)
                           ic = k
                        else
                           ia = k
                           ib = i12(1,k)
                           ic = i
                        end if
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
                        call imager (xcb,ycb,zcb,j)
                        xp = ycb*zab - zcb*yab
                        yp = zcb*xab - xcb*zab
                        zp = xcb*yab - ycb*xab
                        rp = sqrt(xp*xp + yp*yp + zp*zp)
                        rab2 = xab*xab + yab*yab + zab*zab
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        rcb2 = max(0.0001d0,rcb2)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        sine = sqrt(abs(1.0d0-cosine**2))
                        rab = sqrt(rab2)
                        ideal = bl(bndlist(1,ia))
                        ratio = rab / ideal
                        fterm = cosine * ratio
                        deddt = -sine * ratio
                        deddr = cosine / (rab*ideal)
                     end if
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              fterm = 1.0d0
                              rv = radmin4(kt,it)
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     rik = sqrt(rik2)
                     if (p2 .le. expmin2) then
                        e = 0.0d0
                        de = 0.0d0
                     else if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bbuck / rv
                        expterm = abuck * exp(-bbuck/p)
                        fcbuck = fterm * cbuck * p6
                        e = eps * (expterm - fcbuck)
                        de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                     else
                        use_hb = .false.
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        if (use_hb) then
                           deddt = deddt * fgrp
                           deddr = deddr * fgrp
                        end if
                     end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (i .ne. k) then
                        if (k .eq. kv) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kv) = dev(1,kv) - dedx*redkv
                           dev(2,kv) = dev(2,kv) - dedy*redkv
                           dev(3,kv) = dev(3,kv) - dedz*redkv
                        end if
                     end if
c
c     find the chain rule terms for hydrogen bonding components
c
                     if (use_hb) then
                        term = eps * cbuck * p6
                        deddt = deddt * term
                        deddr = deddr * term
                        if (rik2 .gt. cut2) then
                           deddt = deddt * taper
                           deddr = deddr * taper
                        end if
                        terma = deddt / (rab2*max(rp,1.0d-6))
                        termc = -deddt / (rcb2*max(rp,1.0d-6))
                        dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                        dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                        dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                        dedxic = termc * (ycb*zp-zcb*yp)
                        dedyic = termc * (zcb*xp-xcb*zp)
                        dedzic = termc * (xcb*yp-ycb*xp)
                        dedxib = -dedxia - dedxic
                        dedyib = -dedyia - dedyic
                        dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                        dev(1,ia) = dev(1,ia) + dedxia
                        dev(2,ia) = dev(2,ia) + dedyia
                        dev(3,ia) = dev(3,ia) + dedzia
                        dev(1,ib) = dev(1,ib) + dedxib
                        dev(2,ib) = dev(2,ib) + dedyib
                        dev(3,ib) = dev(3,ib) + dedzib
                        dev(1,ic) = dev(1,ic) + dedxic
                        dev(2,ic) = dev(2,ic) + dedyic
                        dev(3,ic) = dev(3,ic) + dedzic
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     if (use_hb) then
                        vxx = xab*dedxia + xcb*dedxic
                        vyx = yab*dedxia + ycb*dedxic
                        vzx = zab*dedxia + zcb*dedxic
                        vyy = yab*dedyia + ycb*dedyic
                        vzy = zab*dedyia + zcb*dedyic
                        vzz = zab*dedzia + zcb*dedzic
                     end if
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                     einter = einter + e
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emm3hb1b  --  MM3 vdw-hbond derivs via lights  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emm3hb1b" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding energy with respect to Cartesian
c     coordinates using the method of lights
c
c
      subroutine emm3hb1b
      use sizes
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use iounit
      use light
      use molcul
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer kgy,kgz
      integer start,stop
      integer, allocatable :: iv14(:)
      real*8 e,de,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rvterm
      real*8 taper,dtaper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge,fcbuck
      real*8 dot,cosine,sine
      real*8 term,terma,termc
      real*8 fterm,ideal,ratio
      real*8 deddr,deddt
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical use_hb
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     special cutoffs for very short and very long range terms
c
      expmin2 = 0.01d0
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (abuck*exp(-bbuck/expcut) - cbuck*(expcut**6))
     &                               / (expcut**12)
c
c     apply any reduction factor to the atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nvdw,xsort,ysort,zsort,unique)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         usei = (use(i) .or. use(iv))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            k = ivdw(kk-((kk-1)/nvdw)*nvdw)
            kv = ired(k)
            prime = (kk .le. nvdw)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            if (proceed)  proceed = (vscale(k) .ne. 0.0d0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_bounds) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k).eq.i .and. prime) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
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
                     if (use_bounds) then
                        if (abs(xcb) .gt. xcell2)
     &                     xcb = xcb - sign(xcell,xcb)
                        if (abs(ycb) .gt. ycell2)
     &                     ycb = ycb - sign(ycell,ycb)
                        if (abs(zcb) .gt. zcell2)
     &                     zcb = zcb - sign(zcell,zcb)
                        if (monoclinic) then
                           xcb = xcb + zcb*beta_cos
                           zcb = zcb * beta_sin
                        else if (triclinic) then
                           xcb = xcb + ycb*gamma_cos + zcb*beta_cos
                           ycb = ycb*gamma_sin + zcb*beta_term
                           zcb = zcb * gamma_term
                        end if
                     end if
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                     deddt = -sine * ratio
                     deddr = cosine / (rab*ideal)
                  end if
                  if (prime)  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                     de = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bbuck / rv
                     expterm = abuck * exp(-bbuck/p)
                     fcbuck = fterm * cbuck * p6
                     e = eps * (expterm - fcbuck)
                     de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (use_hb) then
                        deddt = deddt * fgrp
                        deddr = deddr * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     find the chain rule terms for hydrogen bonding components
c
                  if (use_hb) then
                     term = eps * cbuck * p6
                     deddt = deddt * term
                     deddr = deddr * term
                     if (rik2 .gt. cut2) then
                        deddt = deddt * taper
                        deddr = deddr * taper
                     end if
                     terma = deddt / (rab2*rp)
                     termc = -deddt / (rcb2*rp)
                     dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                     dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                     dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                     dedxic = termc * (ycb*zp-zcb*yp)
                     dedyic = termc * (zcb*xp-xcb*zp)
                     dedzic = termc * (xcb*yp-ycb*xp)
                     dedxib = -dedxia - dedxic
                     dedyib = -dedyia - dedyic
                     dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                     dev(1,ia) = dev(1,ia) + dedxia
                     dev(2,ia) = dev(2,ia) + dedyia
                     dev(3,ia) = dev(3,ia) + dedzia
                     dev(1,ib) = dev(1,ib) + dedxib
                     dev(2,ib) = dev(2,ib) + dedyib
                     dev(3,ib) = dev(3,ib) + dedzib
                     dev(1,ic) = dev(1,ic) + dedxic
                     dev(2,ic) = dev(2,ic) + dedyic
                     dev(3,ic) = dev(3,ic) + dedzic
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  if (use_hb) then
                     vxx = xab*dedxia + xcb*dedxic
                     vyx = yab*dedxia + ycb*dedxic
                     vzx = zab*dedxia + zcb*dedxic
                     vyy = yab*dedyia + ycb*dedyic
                     vzy = zab*dedyia + zcb*dedyic
                     vzz = zab*dedzia + zcb*dedzic
                  end if
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (.not.prime .or. molcule(i).ne.molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine emm3hb1c  --  MM3 vdw-hbond derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "emm3hb1c" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding energy with respect to Cartesian
c     coordinates using a pairwise neighbor list
c
c
      subroutine emm3hb1c
      use sizes
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer, allocatable :: iv14(:)
      real*8 e,de,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rvterm
      real*8 taper,dtaper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 dot,cosine,sine
      real*8 fterm,fcbuck,term
      real*8 deddr,ideal,ratio
      real*8 deddt,terma,termc
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei,use_hb
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     special cutoffs for very short and very long range terms
c
      expmin2 = 0.01d0
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (abuck*exp(-bbuck/expcut) - cbuck*(expcut**6))
     &                               / (expcut**12)
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            k = ivdw(vlst(kk,ii))
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
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
                     call image (xcb,ycb,zcb)
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                     deddt = -sine * ratio
                     deddr = cosine / (rab*ideal)
                  end if
                  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                     de = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bbuck / rv
                     expterm = abuck * exp(-bbuck/p)
                     fcbuck = fterm * cbuck * p6
                     e = eps * (expterm - fcbuck)
                     de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (use_hb) then
                        deddt = deddt * fgrp
                        deddr = deddr * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     find the chain rule terms for hydrogen bonding components
c
                  if (use_hb) then
                     term = eps * cbuck * p6
                     deddt = deddt * term
                     deddr = deddr * term
                     if (rik2 .gt. cut2) then
                        deddt = deddt * taper
                        deddr = deddr * taper
                     end if
                     terma = deddt / (rab2*max(rp,1.0d-6))
                     termc = -deddt / (rcb2*max(rp,1.0d-6))
                     dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                     dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                     dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                     dedxic = termc * (ycb*zp-zcb*yp)
                     dedyic = termc * (zcb*xp-xcb*zp)
                     dedzic = termc * (xcb*yp-ycb*xp)
                     dedxib = -dedxia - dedxic
                     dedyib = -dedyia - dedyic
                     dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                     dev(1,ia) = dev(1,ia) + dedxia
                     dev(2,ia) = dev(2,ia) + dedyia
                     dev(3,ia) = dev(3,ia) + dedzia
                     dev(1,ib) = dev(1,ib) + dedxib
                     dev(2,ib) = dev(2,ib) + dedyib
                     dev(3,ib) = dev(3,ib) + dedzib
                     dev(1,ic) = dev(1,ic) + dedxic
                     dev(2,ic) = dev(2,ic) + dedyic
                     dev(3,ic) = dev(3,ic) + dedzic
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  if (use_hb) then
                     vxx = xab*dedxia + xcb*dedxic
                     vyx = yab*dedxia + ycb*dedxic
                     vzx = zab*dedxia + zcb*dedxic
                     vyy = yab*dedyia + ycb*dedyic
                     vzy = zab*dedyia + zcb*dedyic
                     vzz = zab*dedzia + zcb*dedzic
                  end if
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
