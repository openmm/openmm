c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edipole3  --  dipole-dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edipole3" calculates the dipole-dipole interaction energy;
c     also partitions the energy among the atoms
c
c
      subroutine edipole3
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use dipole
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
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
      logical header,huge
      character*6 mode
c
c
c     zero out the overall dipole interaction energy contribution
c     and set up the constants for the calculation
c
      ned = 0
      ed = 0.0d0
      do i = 1, n
         aed(i) = 0.0d0
      end do
      if (ndipole .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndipole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole-Dipole Interactions :',
     &           //,' Type',15x,'Dipole 1',14x,'Dipole 2',
     &              8x,'Distance',6x,'Energy',/)
      end if
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      mode = 'DIPOLE'
      call switch (mode)
c
c     compute and partition the dipole interaction energy
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
                  ned = ned + 1
                  ed = ed + e
                  aed(i1) = aed(i1) + 0.25d0*e
                  aed(i2) = aed(i2) + 0.25d0*e
                  aed(k1) = aed(k1) + 0.25d0*e
                  aed(k2) = aed(k2) + 0.25d0*e
c
c     increment the total intermolecular energy
c
                  if (molcule(i1) .ne. molcule(k1)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dipole-Dipole',
     &                             ' Interactions :',
     &                          //,' Type',15x,'Dipole 1',14x,
     &                             'Dipole 2',8x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i1,name(i1),i2,name(i2),
     &                                k1,name(k1),k2,name(k2),
     &                                sqrt(r2),e
   30                format (' Dipole',4x,4(i7,'-',a3),f11.4,f12.4)
                  end if
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
                     if (e .ne. 0.0d0)  ned = ned + 1
                     if (i .eq. k) then
                        ed = ed + 0.5d0*e
                        aed(i1) = aed(i1) + 0.25d0*e
                        aed(i2) = aed(i2) + 0.25d0*e
                     else
                        ed = ed + e
                        aed(i1) = aed(i1) + 0.25d0*e
                        aed(i2) = aed(i2) + 0.25d0*e
                        aed(k1) = aed(k1) + 0.25d0*e
                        aed(k2) = aed(k2) + 0.25d0*e
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 10.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Dipole-Dipole',
     &                                ' Interactions :',
     &                             //,' Type',15x,'Dipole 1',14x,
     &                                'Dipole 2',8x,'Distance',
     &                                6x,'Energy',/)
                        end if
                        write (iout,50)  i1,name(i1),i2,name(i2),
     &                                   k1,name(k1),k2,name(k2),
     &                                   sqrt(r2),e
   50                format (' Dipole',4x,4(i7,'-',a3),f11.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
