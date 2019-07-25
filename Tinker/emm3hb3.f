c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb3  --  MM3 vdw & hbond energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb3" calculates the MM3 exp-6 van der Waals and directional
c     charge transfer hydrogen bonding energy, and partitions the energy
c     among the atoms
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
      subroutine emm3hb3
      use sizes
      use analyz
      use atoms
      use energi
      use inform
      use iounit
      use limits
      use vdwpot
      implicit none
      integer i
      real*8 elrc,aelrc
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_lights) then
         call emm3hb3b
      else if (use_vlist) then
         call emm3hb3c
      else
         call emm3hb3a
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
         aelrc = elrc / dble(n)
         do i = 1, n
            aev(i) = aev(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0d0) then
            write (iout,10)  elrc
   10       format (/,' Long Range vdw Correction :',9x,f12.4)
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emm3hb3a  --  double loop MM3 vdw-hb analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emm3hb3" calculates the MM3 exp-6 van der Waals and
c     directional charge transfer hydrogen bonding energy, and
c     partitions the energy among the atoms
c
c
      subroutine emm3hb3a
      use sizes
      use action
      use analyz
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer, allocatable :: iv14(:)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nvdw.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual van der Waals Interactions :',
     &           //,' Type',14x,'Atom Names',20x,'Minimum',
     &              4x,'Actual',6x,'Energy',/)
      end if
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
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
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
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = abuck * exp(-bbuck/p)
                     e = eps * (expterm - fterm*cbuck*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   30                format (' VDW-MM3',3x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
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
                        call imager (xcb,ycb,zcb,j)
                        rab2 = xab*xab + yab*yab + zab*zab
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        rcb2 = max(0.0001d0,rcb2)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        rab = sqrt(rab2)
                        ideal = bl(bndlist(1,ia))
                        fterm = cosine * (rab/ideal)
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
                     if (p2 .le. expmin2) then
                        e = 0.0d0
                     else if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        expterm = abuck * exp(-bbuck/p)
                        e = eps * (expterm - fterm*cbuck*p6)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                                + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component
c
                     if (e .ne. 0.0d0) then
                        nev = nev + 1
                        if (i .eq. k) then
                           ev = ev + 0.5d0*e
                           aev(i) = aev(i) + 0.5d0*e
                        else
                           ev = ev + e
                           aev(i) = aev(i) + 0.5d0*e
                           aev(k) = aev(k) + 0.5d0*e
                        end if
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                     huge = (e .gt. 10.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual van der Waals',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                20x,'Minimum',4x,'Actual',
     &                                6x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),
     &                                   rv,sqrt(rik2),e
   50                   format (' VDW-MM3',3x,2(i7,'-',a3),3x,
     &                             '(XTAL)',4x,2f10.4,f12.4)
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
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emm3hb3b  --  MM3 vdw-hbond analysis via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emm3hb3b" calculates the MM3 exp-6 van der Waals and
c     directional charge transfer hydrogen bonding energy using
c     the method of lights
c
c
      subroutine emm3hb3b
      use sizes
      use action
      use analyz
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use light
      use molcul
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer kgy,kgz
      integer start,stop
      integer ikmin,ikmax
      integer, allocatable :: iv14(:)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical header,huge
      character*6 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nvdw.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual van der Waals Interactions :',
     &           //,' Type',14x,'Atom Names',20x,'Minimum',
     &              4x,'Actual',6x,'Energy',/)
      end if
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
   20    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 60
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 60
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 60
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 60
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
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k).eq.i .and. prime) then
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
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  if (prime)  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = abuck * exp(-bbuck/p)
                     e = eps * (expterm - fterm*cbuck*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (.not.prime .or. molcule(i).ne.molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     ikmin = min(i,k)
                     ikmax = max(i,k)
                     if (prime) then
                        write (iout,40)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),rv,sqrt(rik2),e
   40                   format (' VDW-MM3',3x,2(i7,'-',a3),
     &                             13x,2f10.4,f12.4)
                     else
                        write (iout,50)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),rv,sqrt(rik2),e
   50                   format (' VDW-MM3',3x,2(i7,'-',a3),3x,
     &                             '(XTAL)',4x,2f10.4,f12.4)
                     end if
                  end if
               end if
            end if
   60       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 20
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emm3hb3c  --  MM3 vdw-hbond analysis via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emm3hb3c" calculates the MM3 exp-6 van der Waals and
c     directional charge transfer hydrogen bonding energy using
c     a pairwise neighbor list
c
c
      subroutine emm3hb3c
      use sizes
      use action
      use analyz
      use atmlst
      use atomid
      use atoms
      use bndstr
      use bound
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer ia,ib,ic
      integer, allocatable :: iv14(:)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmin2
      real*8 expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ev = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
      end do
      if (nvdw .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nvdw.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual van der Waals Interactions :',
     &           //,' Type',14x,'Atom Names',20x,'Minimum',
     &              4x,'Actual',6x,'Energy',/)
      end if
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
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
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
                     rab2 = xab*xab + yab*yab + zab*zab
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     rcb2 = max(0.0001d0,rcb2)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     rab = sqrt(rab2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  eps = eps * vscale(k)
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expmin2) then
                     e = 0.0d0
                  else if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = abuck * exp(-bbuck/p)
                     e = eps * (expterm - fterm*cbuck*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (e .ne. 0.0d0) then
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   30                format (' VDW-MM3',3x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
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
