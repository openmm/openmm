c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echgdpl  --  charge-dipole potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echgdpl" calculates the charge-dipole interaction energy
c
c
      subroutine echgdpl
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use dipole
      use energi
      use group
      use shunt
      use units
      use usage
      implicit none
      integer i,j,k
      integer i1,k1,k2
      integer, allocatable :: skip(:)
      real*8 e,rk2,rkr3,dotk
      real*8 taper,fgrp
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      logical proceed
      character*6 mode
c
c
c     zero out the overall charge-dipole interaction energy
c
      ecd = 0.0d0
      if (nion.eq.0 .or. ndipole.eq.0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (skip(n))
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      mode = 'CHGDPL'
      call switch (mode)
c
c     get the total energy by looping over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = x(k1) + xk*sdpl(k) - xi
               yr = y(k1) + yk*sdpl(k) - yi
               zr = z(k1) + zk*sdpl(k) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * bdpl(k)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                  ecd = ecd + e
               end if
            end if
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
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = x(k1) + xk*sdpl(k) - xi
                  yr = y(k1) + yk*sdpl(k) - yi
                  zr = z(k1) + zk*sdpl(k) - zi
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * bdpl(k)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (skip(k1).eq.i1 .or. skip(k2).ne.i1)
     &                        fik = 0.0d0
                        end if
                     end if
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                     ecd = ecd + e
                  end if
               end do
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (skip)
      return
      end
