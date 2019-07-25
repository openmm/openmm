c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emm3hb2  --  atomwise MM3 vdw & hbond Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emm3hb2" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding second derivatives for a single
c     atom at a time
c
c     note this version only partially incorporates the directional
c     hydrogen bonding term into the Hessian calculation
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
      subroutine emm3hb2 (iatom,xred,yred,zred)
      use sizes
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use cell
      use chgpot
      use couple
      use group
      use hessn
      use shunt
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer iatom,jcell
      integer ia,ib,ic
      integer nlist,list(5)
      integer, allocatable :: iv14(:)
      real*8 e,de,d2e
      real*8 eps,rv,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk
      real*8 redikv,redivkv
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper,d2taper
      real*8 d2edx,d2edy,d2edz
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 rvterm,rvterm2
      real*8 dot,cosine,sine
      real*8 fterm,fcbuck
      real*8 ideal,ratio
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 term(3,3)
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
      real*8, allocatable :: vscale(:)
      logical proceed
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
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
c     check to see if the atom of interest is a vdw site
c
      nlist = 0
      do k = 1, nvdw
         if (ivdw(k) .eq. iatom) then
            nlist = nlist + 1
            list(nlist) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     determine the atoms involved via reduction factors
c
      nlist = 1
      list(nlist) = iatom
      do k = 1, n12(iatom)
         i = i12(k,iatom)
         if (ired(i) .eq. iatom) then
            nlist = nlist + 1
            list(nlist) = i
         end if
      end do
c
c     find van der Waals Hessian elements for involved atoms
c
      do ii = 1, nlist
         i = list(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
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
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (k .ne. i)
c
c     compute the Hessian elements for this interaction
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
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
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
                  end if
                  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                     de = 0.0d0
                     d2e = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bbuck / rv
                     rvterm2 = rvterm * rvterm
                     expterm = abuck * exp(-bbuck/p)
                     fcbuck = fterm * cbuck * p6
                     e = eps * (expterm - fcbuck)
                     de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                     d2e = eps * (rvterm2*expterm-42.0d0*fcbuck/rik2)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                     d2e = 156.0d0 * e / rik2
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
                     d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                            + 6.0d0*c3*rik + 2.0d0*c2
                     d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                     de = e*dtaper + de*taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                  de = de / rik
                  d2e = (d2e-de) / rik2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  if (i .eq. iatom) then
                     if (i.eq.iv .and. k.eq.kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)
                           hessy(j,k) = hessy(j,k) - term(2,j)
                           hessz(j,k) = hessz(j,k) - term(3,j)
                        end do
                     else if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redi
                           hessy(j,k) = hessy(j,k) - term(2,j)*redi
                           hessz(j,k) = hessz(j,k) - term(3,j)*redi
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                        end do
                     else if (i .eq. iv) then
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)*redk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redk
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redkv
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redik = redi * redk
                        redikv = redi * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redik
                           hessy(j,k) = hessy(j,k) - term(2,j)*redik
                           hessz(j,k) = hessz(j,k) - term(3,j)*redik
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redikv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redikv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redikv
                        end do
                     end if
                  else if (iv .eq. iatom) then
                     if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                           hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                           hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redivk = rediv * redk
                        redivkv = rediv * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redivkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redivkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redivkv
                        end do
                     end if
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
      do ii = 1, nlist
         i = list(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
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
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the Hessian elements for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               do jcell = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call imager (xr,yr,zr,jcell)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
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
                        call imager (xcb,ycb,zcb,jcell)
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
                        d2e = 0.0d0
                     else if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bbuck / rv
                        rvterm2 = rvterm * rvterm
                        expterm = abuck * exp(-bbuck/p)
                        fcbuck = fterm * cbuck * p6
                        e = eps * (expterm - fcbuck)
                        de = eps * (rvterm*expterm+6.0d0*fcbuck/rik)
                        d2e = eps * (rvterm2*expterm-42.0d0*fcbuck/rik2)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                        d2e = 156.0d0 * e / rik2
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
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                             + 6.0d0*c3*rik + 2.0d0*c2
                        d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                        de = e*dtaper + de*taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        de = de * fgrp
                        d2e = d2e * fgrp
                     end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                     de = de / rik
                     d2e = (d2e-de) / rik2
                     d2edx = d2e * xr
                     d2edy = d2e * yr
                     d2edz = d2e * zr
                     term(1,1) = d2edx*xr + de
                     term(1,2) = d2edx*yr
                     term(1,3) = d2edx*zr
                     term(2,1) = term(1,2)
                     term(2,2) = d2edy*yr + de
                     term(2,3) = d2edy*zr
                     term(3,1) = term(1,3)
                     term(3,2) = term(2,3)
                     term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                     if (i .eq. iatom) then
                        if (i.eq.iv .and. k.eq.kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)
                              hessy(j,k) = hessy(j,k) - term(2,j)
                              hessz(j,k) = hessz(j,k) - term(3,j)
                           end do
                        else if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redi
                              hessy(j,k) = hessy(j,k) - term(2,j)*redi
                              hessz(j,k) = hessz(j,k) - term(3,j)*redi
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                           end do
                        else if (i .eq. iv) then
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)*redk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redk
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redkv
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redik = redi * redk
                           redikv = redi * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redik
                              hessy(j,k) = hessy(j,k) - term(2,j)*redik
                              hessz(j,k) = hessz(j,k) - term(3,j)*redik
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redikv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redikv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redikv
                           end do
                        end if
                     else if (iv .eq. iatom) then
                        if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                              hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                              hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redivk = rediv * redk
                           redivkv = rediv * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redivkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redivkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redivkv
                           end do
                        end if
                     end if
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
      deallocate (vscale)
      return
      end
