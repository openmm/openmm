c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine edipole  --  dipole-dipole potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "edipole" calculates the dipole-dipole interaction energy
c
c
      subroutine edipole
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use dipole
      use energi
      use group
      use shunt
      use units
      use usage
      implicit none
      integer i,j,k
      integer i1,i2,k1,k2
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xq,yq,zq
      real*8 xr,yr,zr
      real*8 f,fi,fik
      real*8 taper,fgrp
      real*8 e,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 r,r2,r3,r4,r5
      logical proceed
      character*6 mode
c
c
c     zero out the overall dipole interaction energy
c     and set up the constants for the calculation
c
      ed = 0.0d0
      if (ndipole .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      mode = 'DIPOLE'
      call switch (mode)
c
c     calculate the pairwise dipole interaction energy term
c
      do i = 1, ndipole-1
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
c
c     decide whether to compute the current interaction
c
         do k = i+1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               if (use_polymer)  call imager (xk,yk,zk,-1)
               xr = xq - x(k1) - xk*sdpl(k)
               yr = yq - y(k1) - yk*sdpl(k)
               zr = zq - z(k1) - zk*sdpl(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
                  e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                  ed = ed + e
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
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
c
c     decide whether to compute the current interaction
c
         do k = i, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  if (use_polymer)  call imager (xk,yk,zk,-1)
                  xr = xq - x(k1) - xk*sdpl(k)
                  yr = yq - y(k1) - yk*sdpl(k)
                  zr = zq - z(k1) - zk*sdpl(k)
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rirkr3 = sqrt(ri2*rk2*r2) * r2
                     dotp = xi*xk + yi*yk + zi*zk
                     doti = xi*xr + yi*yr + zi*zr
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fi * bdpl(k)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (k1.eq.i1 .or. k1.eq.i2 .or.
     &                         k2.eq.i1 .or. k2.eq.i2)  fik = 0.0d0
                        end if
                     end if
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ed = ed + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
