c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine egauss3  --  Gaussian vdw energy & analysis  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egauss3" calculates the Gaussian expansion van der Waals
c     interaction energy and partitions the energy among the atoms
c
c
      subroutine egauss3
      use sizes
      use limits
      use warp
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call egauss3d
      else if (use_vlist) then
         call egauss3c
      else if (use_lights) then
         call egauss3b
      else
         call egauss3a
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine egauss3a  --  double loop Gaussian vdw analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "egauss3a" calculates the Gaussian expansion van der Waals
c     energy and partitions the energy among the atoms using a
c     pairwise double loop
c
c
      subroutine egauss3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
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
      integer i,j,k,m
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 rad,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expterm
      real*8 a(maxgauss)
      real*8 b(maxgauss)
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
      expcut = -50.0d0
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
                  rad2 = radmin(kt,it)**2
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rad2 = radmin4(kt,it)**2
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                  end do
                  e = 0.0d0
                  do j = 1, ngauss
                     expterm = -b(j) * rik2
                     if (expterm .gt. expcut)
     &                  e = e + a(j)*exp(expterm)
                  end do
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
c     increment the total van der Waals energy components
c
                  nev = nev + 1
                  aev(i) = aev(i) + 0.5d0*e
                  aev(k) = aev(k) + 0.5d0*e
                  ev = ev + e
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
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     rad = sqrt(rad2)
                     write (iout,30)  i,name(i),k,name(k),
     &                                rad,sqrt(rik2),e
   30                format (' VDW-Gauss',1x,2(i7,'-',a3),
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
               do m = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call imager (xr,yr,zr,m)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     rad2 = radmin(kt,it)**2
                     eps = epsilon(kt,it)
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              rad2 = radmin4(kt,it)**2
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
                     do j = 1, ngauss
                        a(j) = igauss(1,j) * eps
                        b(j) = igauss(2,j) / rad2
                     end do
                     e = 0.0d0
                     do j = 1, ngauss
                        expterm = -b(j) * rik2
                        if (expterm .gt. expcut)
     &                     e = e + a(j)*exp(expterm)
                     end do
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
                        rad = sqrt(rad2)
                        write (iout,50)  i,name(i),k,name(k),
     &                                   rad,sqrt(rik2),e
   50                   format (' VDW-Gauss',1x,2(i7,'-',a3),
     &                             3x,'(XTAL)',4x,2f10.4,f12.4)
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egauss3b  --  Gaussian vdw analysis via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss3b" calculates the Gaussian expansion van der Waals
c     energy and partitions the energy among the atoms using the
c     method of lights
c
c
      subroutine egauss3b
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
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
      integer i,j,k,m
      integer ii,iv,it
      integer kk,kv,kt
      integer kgy,kgz
      integer start,stop
      integer ikmin,ikmax
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 rad,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expterm
      real*8 a(maxgauss)
      real*8 b(maxgauss)
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
c     set cutoff distances and switching function coefficients
c
      mode = 'VDW'
      call switch (mode)
      expcut = -50.0d0
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
c     loop over all atoms computing the interactions
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
         do m = start, stop
            kk = locx(m)
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
               xr = xi - xsort(m)
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
                  rad2 = radmin(kt,it)**2
                  eps = epsilon(kt,it)
                  if (prime) then
                     if (iv14(k) .eq. i) then
                        rad2 = radmin4(kt,it)**2
                        eps = epsilon4(kt,it)
                     end if
                     eps = eps * vscale(k)
                  end if
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                  end do
                  e = 0.0d0
                  do j = 1, ngauss
                     expterm = -b(j) * rik2
                     if (expterm .gt. expcut)
     &                  e = e + a(j)*exp(expterm)
                  end do
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
                     rad = sqrt(rad2)
                     ikmin = min(i,k)
                     ikmax = max(i,k)
                     if (prime) then
                        write (iout,40)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),rad,sqrt(rik2),e
   40                   format (' VDW-Gauss',1x,2(i7,'-',a3),
     &                             13x,2f10.4,f12.4)
                     else
                        write (iout,50)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),rad,sqrt(rik2),e
   50                   format (' VDW-Gauss',1x,2(i7,'-',a3),
     &                             3x,'(XTAL)',4x,2f10.4,f12.4)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egauss3c  --  Gaussian vdw analysis via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "egauss3c" calculates the Gaussian expansion van der Waals
c     energy and partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine egauss3c
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
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
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 rad,rad2,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expterm
      real*8 a(maxgauss)
      real*8 b(maxgauss)
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
c     set cutoff distances and switching function coefficients
c
      mode = 'VDW'
      call switch (mode)
      expcut = -50.0d0
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
                  rad2 = radmin(kt,it)**2
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rad2 = radmin4(kt,it)**2
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                  end do
                  e = 0.0d0
                  do j = 1, ngauss
                     expterm = -b(j) * rik2
                     if (expterm .gt. expcut)
     &                  e = e + a(j)*exp(expterm)
                  end do
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
c     increment the total van der Waals energy components
c
                  nev = nev + 1
                  aev(i) = aev(i) + 0.5d0*e
                  aev(k) = aev(k) + 0.5d0*e
                  ev = ev + e
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
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             20x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     rad = sqrt(rad2)
                     write (iout,30)  i,name(i),k,name(k),
     &                                rad,sqrt(rik2),e
   30                format (' VDW-Gauss',1x,2(i7,'-',a3),
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egauss3d  --  Gaussian analysis for smoothing  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egauss3d" calculates the Gaussian expansion van der Waals
c     interaction energy and partitions the energy among the atoms
c     for use with potential energy smoothing
c
c
      subroutine egauss3d
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use usage
      use vdw
      use vdwpot
      use warp
      implicit none
      integer i,j,k,ii,kk
      integer iv,kv,it,kt
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 rad,rad2,fgrp
      real*8 rik,rik2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,expcut,broot
      real*8 expterm,expterm2
      real*8 width,wterm
      real*8 t1,t2,term
      real*8 a(maxgauss)
      real*8 b(maxgauss)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical header,huge
      external erf
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
c     set the extent of smoothing to be performed
c
      expcut = -50.0d0
      width = 0.0d0
      if (use_dem) then
         width = 4.0d0 * diffv * deform
      else if (use_gda) then
         wterm = (2.0d0/3.0d0) * diffv
      else if (use_tophat) then
         width = max(diffv*deform,0.0001d0)
      end if
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
         do j = ii+1, nvdw
            vscale(ivdw(j)) = 1.0d0
         end do
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
            if (proceed)  proceed = (vscale(k) .ne. 0.0d0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               rad2 = radmin(kt,it)**2
               eps = epsilon(kt,it)
               if (iv14(k) .eq. i) then
                  rad2 = radmin4(kt,it)**2
                  eps = epsilon4(kt,it)
               end if
               eps = eps * vscale(k)
               do j = 1, ngauss
                  a(j) = igauss(1,j) * eps
                  b(j) = igauss(2,j) / rad2
               end do
               e = 0.0d0
c
c     transform the potential function via smoothing
c
               if (use_tophat) then
                  rik = sqrt(rik2)
                  do j = 1, ngauss
                     broot = sqrt(b(j))
                     expterm = -b(j) * (rik+width)**2
                     if (expterm .gt. expcut) then
                        expterm = exp(expterm)
                     else
                        expterm = 0.0d0
                     end if
                     expterm2 = -b(j) * (width-rik)**2
                     if (expterm2 .gt. expcut) then
                        expterm2 = exp(expterm2)
                     else
                        expterm2 = 0.0d0
                     end if
                     term = broot * (expterm-expterm2)
                     term = term + sqrtpi*b(j)*rik
     &                         * (erf(broot*(rik+width))
     &                           +erf(broot*(width-rik)))
                     e = e + term*a(j)/(b(j)*b(j)*broot)
                  end do
                  e = e * 3.0d0/(8.0d0*rik*width**3)
               else
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut)
     &                  e = e + (a(j)/t2)*exp(expterm)
                  end do
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total van der Waals energy components
c
               nev = nev + 1
               aev(i) = aev(i) + 0.5d0*e
               aev(k) = aev(k) + 0.5d0*e
               ev = ev + e
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
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual van der Waals',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          20x,'Minimum',4x,'Actual',
     &                          6x,'Energy',/)
                  end if
                  rad = sqrt(rad2)
                  write (iout,30)  i,name(i),k,name(k),
     &                             rad,sqrt(rik2),e
   30             format (' VDW-Gauss',1x,2(i7,'-',a3),
     &                       13x,2f10.4,f12.4)
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
