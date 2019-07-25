c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3  --  buffered 14-7 vdw energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms
c
c
      subroutine ehal3
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
         call ehal3b
      else if (use_vlist) then
         call ehal3c
      else
         call ehal3a
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal3a  --  double loop buffered 14-7 analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal3a" calculates the buffered 14-7 van der Waals energy
c     and partitions the energy among the atoms using a pairwise
c     double loop
c
c
      subroutine ehal3a
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
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical muti,mutk
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
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
   30                format (' VDW-Hal',3x,2(i7,'-',a3),
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
         muti = mut(i)
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
            mutk = mut(k)
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
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              rv = radmin4(kt,it)
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
c
c     get the interaction energy, via soft core if necessary
c
                     if ((muti .and. .not.mutk) .or.
     &                   (mutk .and. .not.muti)) then
                        rho = rik / rv
                        eps = eps * vlambda**scexp
                        scal = scalpha * (1.0d0-vlambda)**2
                        t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                        t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                        e = eps * t1 * (t2-2.0d0)
                     else
                        rv7 = rv**7
                        rik7 = rik**7
                        rho = rik7 + ghal*rv7
                        tau = (dhal+1.0d0) / (rik + dhal*rv)
                        e = eps * rv7 * tau**7
     &                         * ((ghal+1.0d0)*rv7/rho-2.0d0)
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
   50                   format (' VDW-Hal',3x,2(i7,'-',a3),3x,
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal3b  --  buffered 14-7 analysis via lights  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal3b" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using the
c     method of lights
c
c
      subroutine ehal3b
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
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer kgy,kgz
      integer start,stop
      integer ikmin,ikmax
      integer, allocatable :: iv14(:)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical muti,mutk
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (prime) then
                     if (iv14(k) .eq. i) then
                        rv = radmin4(kt,it)
                        eps = epsilon4(kt,it)
                     end if
                     eps = eps * vscale(k)
                  end if
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
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
   40                   format (' VDW-Hal',3x,2(i7,'-',a3),
     &                             13x,2f10.4,f12.4)
                     else
                        write (iout,50)  ikmax,name(ikmin),ikmax,
     &                                   name(ikmax),rv,sqrt(rik2),e
   50                   format (' VDW-Hal',3x,2(i7,'-',a3),3x,
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ehal3c  --  buffered 14-7 analysis via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ehal3c" calculates the buffered 14-7 van der Waals energy
c     and also partitions the energy among the atoms using a
c     pairwise neighbor list
c
c
      subroutine ehal3c
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
      use mutant
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
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical muti,mutk
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
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
!$OMP& jvdw,xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
!$OMP& use_group,off2,radmin,epsilon,radmin4,epsilon4,ghal,dhal,
!$OMP& vlambda,scalpha,scexp,mut,cut2,c0,c1,c2,c3,c4,c5,
!$OMP& molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(vscale,iv14) shared(ev,einter,nev,aev)
!$OMP DO reduction(+:ev,einter,nev,aev) schedule(guided)
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
c
c     get the interaction energy, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
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
   30                format (' VDW-Hal',3x,2(i7,'-',a3),
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
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
