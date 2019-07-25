c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge  --  charge-charge potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge" calculates the charge-charge interaction energy
c
c
      subroutine echarge
      use sizes
      use limits
      use warp
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call echarge0g
      else if (use_ewald) then
         if (use_clist) then
            call echarge0f
         else if (use_lights) then
            call echarge0e
         else
            call echarge0d
         end if
      else if (use_clist) then
         call echarge0c
      else if (use_lights) then
         call echarge0b
      else
         call echarge0a
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge0a  --  charge energy via double loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge0a" calculates the charge-charge interaction energy
c     using a pairwise double loop
c
c
      subroutine echarge0a
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use group
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer in,kn
      integer ic,kc
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     calculate the charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         usei = (use(i) .or. use(ic))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
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
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         usei = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call imager (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     rb = r + ebuffer
                     fik = fi * pchg(kk)
                     if (use_polymer) then
                        if (r2 .le. polycut2)  fik = fik * cscale(kn)
                     end if
                     e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        e = e*taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0b  --  method of lights charge energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0b" calculates the charge-charge interaction energy
c     using the method of lights
c
c
      subroutine echarge0b
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use group
      use iounit
      use light
      use shunt
      use usage
      implicit none
      integer i,j,k,ii,kk
      integer in,ic,kn,kc
      integer kgy,kgz,kmap
      integer start,stop
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      character*6 mode
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = kion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(ii)
         usei = (use(i) .or. use(ic))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
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
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(kmap)
            kc = kion(kmap)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_bounds) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kmap)
                  if (prime)  fik = fik * cscale(kn)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(kmap)
                  end if
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
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
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge0c  --  charge energy via neighbor list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge0c" calculates the charge-charge interaction energy
c     using a pairwise neighbor list
c
c
      subroutine echarge0c
      use sizes
      use atoms
      use bound
      use charge
      use chgpot
      use couple
      use energi
      use group
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer in,kn,ic,kc
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,kion,use,
!$OMP& x,y,z,f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,
!$OMP& c2scale,c3scale,c4scale,c5scale,use_group,use_bounds,off,off2,
!$OMP& cut,cut2,c0,c1,c2,c3,c4,c5,f0,f1,f2,f3,f4,f5,f6,f7,ebuffer)
!$OMP& firstprivate(cscale) shared (ec)
!$OMP DO reduction(+:ec) schedule(guided)
c
c     calculate the charge interaction energy term
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         usei = (use(i) .or. use(ic))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
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
      deallocate (cscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge0d  --  double loop Ewald charge energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge0d" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation
c
c
      subroutine echarge0d
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use math
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,fs,fgrp
      real*8 f,fi,fik
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      character*6 mode
      external erfc
c
c
c     zero out the Ewald charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  scaleterm = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate real space portion involving other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rb = r + ebuffer
                     fik = fi * pchg(kk)
                     rew = aewald * r
                     erfterm = erfc (rew)
                     scale = 1.0d0
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * cscale(kn)
                        end if
                     end if
                     scaleterm = scale - 1.0d0
                     e = (fik/rb) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge0e  --  Ewald charge energy via lights  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge0e" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation and the method of lights
c
c
      subroutine echarge0e
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use light
      use math
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,in,kn
      integer kgy,kgz,kmap
      integer start,stop
      real*8 e,fs,fgrp
      real*8 f,fi,fik
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      character*6 mode
      external erfc
c
c
c     zero out the Ewald charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip
c
c     compute the real space portion of the Ewald summation;
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = iion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         fi = f * pchg(ii)
         usei = (use(i))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
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
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(kmap)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
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
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = 1.0d0
                  if (prime)  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  fik = fi * pchg(kmap)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(kmap)
                  end if
                  scaleterm = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
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
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge0f  --  Ewald charge energy via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge0f" calculates the charge-charge interaction energy
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge0f
      use sizes
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use math
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer in,kn
      real*8 e,fs,fgrp
      real*8 f,fi,fik
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale,scaleterm
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      character*6 mode
      external erfc
c
c
c     zero out the Ewald charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / sqrtpi
      do ii = 1, nion
         e = fs * pchg(ii)**2
         ec = ec + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,use,x,y,z,
!$OMP& f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,c2scale,
!$OMP& c3scale,c4scale,c5scale,use_group,off2,aewald,ebuffer)
!$OMP& firstprivate(cscale) shared (ec)
!$OMP DO reduction(+:ec) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  scaleterm = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scaleterm)
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
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
      deallocate (cscale)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine echarge0g  --  charge energy for smoothing  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "echarge0g" calculates the charge-charge interaction energy
c     for use with potential smoothing methods
c
c
      subroutine echarge0g
      use sizes
      use atoms
      use charge
      use chgpot
      use couple
      use energi
      use group
      use usage
      use warp
      implicit none
      integer i,j,k
      integer ii,kk
      integer in,kn
      real*8 e,fgrp
      real*8 r,r2,rb,rb2
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,wterm,width
      real*8 width2,width3
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      external erf
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set the energy units conversion factor
c
      f = electric / dielec
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     calculate the charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         usei = (use(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rb = r + ebuffer
               fik = fi * pchg(kk) * cscale(kn)
               e = fik / rb
c
c     transform the potential function via smoothing
c
               if (use_dem) then
                  if (width .gt. 0.0d0) then
                     e = e * erf(width*rb)
                  end if
               else if (use_gda) then
                  width = m2(i) + m2(k)
                  if (width .gt. 0.0d0) then
                     width = wterm / sqrt(width)
                     e = e * erf(width*rb)
                  end if
               else if (use_tophat) then
                  if (width .gt. rb) then
                     rb2 = rb * rb
                     e = fik * (3.0d0*width2-rb2) / (2.0d0*width3)
                  end if
               else if (use_stophat) then
                  e = fik / (rb+width)
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
               ec = ec + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ecrecip  --  PME reciprocal space charge energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ecrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to partial charges
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip
      use sizes
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use pme
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      real*8 e,f,denom
      real*8 term,expterm
      real*8 pterm,volterm
      real*8 hsq,struc2
      real*8 h1,h2,h3
      real*8 r1,r2,r3
c
c
c     zero out the particle mesh Ewald grid values
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     get B-spline coefficients and put charges onto grid
c
      call bspline_fill
      call table_fill
      call grid_pchg
c
c     perform the 3-D FFT forward transformation
c
      call fftfront
c
c     use scalar sum to get the reciprocal space energy
c
      f = 0.5d0 * electric / dielec
      npoint = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, npoint-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = f * expterm * struc2
            ec = ec + e
         end if
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         ec = ec + e
      end if
      return
      end
