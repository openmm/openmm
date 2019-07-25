c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdrest  --  stop system translation & rotation  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdrest" finds and removes any translational or rotational
c     kinetic energy of the overall system center of mass
c
c
      subroutine mdrest (istep)
      use sizes
      use atomid
      use atoms
      use bound
      use group
      use inform
      use iounit
      use mdstuf
      use moldyn
      use rgddyn
      use units
      implicit none
      integer i,j,k,istep
      real*8 etrans,erot
      real*8 weigh,totmass,eps
      real*8 xx,yy,zz,xy,xz,yz
      real*8 xtot,ytot,ztot
      real*8 xdel,ydel,zdel
      real*8 mang(3),vang(3)
      real*8 vtot(3),tensor(3,3)
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
c
c
c     check steps between center of mass motion removal
c
      if (.not.dorest)  return
      if (mod(istep,irest) .ne. 0)  return
c
c     zero out the total mass and overall linear velocity
c
      totmass = 0.0d0
      do j = 1, 3
         vtot(j) = 0.0d0
      end do
c
c     compute linear velocity of the system center of mass
c
      if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            weigh = grpmass(i)
            totmass = totmass + weigh
            do j = 1, 3
               vtot(j) = vtot(j) + vcm(j,i)*weigh
            end do
         end do
      else
         do i = 1, n
            weigh = mass(i)
            totmass = totmass + weigh
            do j = 1, 3
               vtot(j) = vtot(j) + v(j,i)*weigh
            end do
         end do
      end if
c
c     compute translational kinetic energy of overall system
c
      etrans = 0.0d0
      do j = 1, 3
         vtot(j) = vtot(j) / totmass
         etrans = etrans + vtot(j)**2
      end do
      etrans = 0.5d0 * etrans * totmass / convert
c
c     perform dynamic allocation of some local arrays
c
      if (.not.use_bounds .and. integrate.eq.'RIGIDBODY') then
         allocate (xcm(ngrp))
         allocate (ycm(ngrp))
         allocate (zcm(ngrp))
      end if
c
c     find the center of mass coordinates of the overall system
c
      if (.not. use_bounds) then
         xtot = 0.0d0
         ytot = 0.0d0
         ztot = 0.0d0
         if (integrate .eq. 'RIGIDBODY') then
            do i = 1, ngrp
               xcm(i) = 0.0d0
               ycm(i) = 0.0d0
               zcm(i) = 0.0d0
               do j = igrp(1,i), igrp(2,i)
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm(i) = xcm(i) + x(k)*weigh
                  ycm(i) = ycm(i) + y(k)*weigh
                  zcm(i) = zcm(i) + z(k)*weigh
               end do
               xtot = xtot + xcm(i)
               ytot = ytot + ycm(i)
               ztot = ztot + zcm(i)
               weigh = max(1.0d0,grpmass(i))
               xcm(i) = xcm(i) / weigh
               ycm(i) = ycm(i) / weigh
               zcm(i) = zcm(i) / weigh
            end do
         else
            do i = 1, n
               weigh = mass(i)
               xtot = xtot + x(i)*weigh
               ytot = ytot + y(i)*weigh
               ztot = ztot + z(i)*weigh
            end do
         end if
         xtot = xtot / totmass
         ytot = ytot / totmass
         ztot = ztot / totmass
c
c     compute the angular momentum of the overall system
c
         do j = 1, 3
            mang(j) = 0.0d0
         end do
         if (integrate .eq. 'RIGIDBODY') then
            do i = 1, ngrp
               weigh = grpmass(i)
               mang(1) = mang(1) + (ycm(i)*vcm(3,i)
     &                             -zcm(i)*vcm(2,i))*weigh
               mang(2) = mang(2) + (zcm(i)*vcm(1,i)
     &                             -xcm(i)*vcm(3,i))*weigh
               mang(3) = mang(3) + (xcm(i)*vcm(2,i)
     &                             -ycm(i)*vcm(1,i))*weigh
            end do
         else
            do i = 1, n
               weigh = mass(i)
               mang(1) = mang(1) + (y(i)*v(3,i)-z(i)*v(2,i))*weigh
               mang(2) = mang(2) + (z(i)*v(1,i)-x(i)*v(3,i))*weigh
               mang(3) = mang(3) + (x(i)*v(2,i)-y(i)*v(1,i))*weigh
            end do
         end if
         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
c
c     calculate the moment of inertia tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         if (integrate .eq. 'RIGIDBODY') then
            do i = 1, ngrp
               weigh = grpmass(i)
               xdel = xcm(i) - xtot
               ydel = ycm(i) - ytot
               zdel = zcm(i) - ztot
               xx = xx + xdel*xdel*weigh
               xy = xy + xdel*ydel*weigh
               xz = xz + xdel*zdel*weigh
               yy = yy + ydel*ydel*weigh
               yz = yz + ydel*zdel*weigh
               zz = zz + zdel*zdel*weigh
            end do
         else
            do i = 1, n
               weigh = mass(i)
               xdel = x(i) - xtot
               ydel = y(i) - ytot
               zdel = z(i) - ztot
               xx = xx + xdel*xdel*weigh
               xy = xy + xdel*ydel*weigh
               xz = xz + xdel*zdel*weigh
               yy = yy + ydel*ydel*weigh
               yz = yz + ydel*zdel*weigh
               zz = zz + zdel*zdel*weigh
            end do
         end if
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
c
c     fix to avoid singularity for one- or two-body systems
c
         if (integrate .eq. 'RIGIDBODY') then
            if (ngrp .le. 2) then
               eps = 0.000001d0
               tensor(1,1) = tensor(1,1) + eps
               tensor(2,2) = tensor(2,2) + eps
               tensor(3,3) = tensor(3,3) + eps
            end if
         else
            if (n .le. 2) then
               eps = 0.000001d0
               tensor(1,1) = tensor(1,1) + eps
               tensor(2,2) = tensor(2,2) + eps
               tensor(3,3) = tensor(3,3) + eps
            end if
         end if
c
c     diagonalize the moment of inertia tensor
c
         call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0d0
         do i = 1, 3
            vang(i) = 0.0d0
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5d0 * erot / convert
      end if
c
c     eliminate any translation of the overall system
c
      if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            do j = 1, 3
               vcm(j,i) = vcm(j,i) - vtot(j)
            end do
         end do
      else
         do i = 1, n
            do j = 1, 3
               v(j,i) = v(j,i) - vtot(j)
            end do
         end do
      end if
c
c     print the translational velocity of the overall system
c
      if (debug) then
         write (iout,10)  (vtot(i),i=1,3),etrans
   10    format (' System Linear Velocity :  ',3d12.2,
     &           /,' Translational Kinetic Energy :',10x,f12.4,
     &              ' Kcal/mole')
      end if
c
c     eliminate any rotation about the system center of mass
c
      if (.not. use_bounds) then
         if (integrate .eq. 'RIGIDBODY') then
            do i = 1, ngrp
               xdel = xcm(i) - xtot
               ydel = ycm(i) - ytot
               zdel = zcm(i) - ztot
               vcm(1,i) = vcm(1,i) - vang(2)*zdel + vang(3)*ydel
               vcm(2,i) = vcm(2,i) - vang(3)*xdel + vang(1)*zdel
               vcm(3,i) = vcm(3,i) - vang(1)*ydel + vang(2)*xdel
            end do
         else
            do i = 1, n
               xdel = x(i) - xtot
               ydel = y(i) - ytot
               zdel = z(i) - ztot
               v(1,i) = v(1,i) - vang(2)*zdel + vang(3)*ydel
               v(2,i) = v(2,i) - vang(3)*xdel + vang(1)*zdel
               v(3,i) = v(3,i) - vang(1)*ydel + vang(2)*xdel
            end do
         end if
c
c     print the angular velocity of the overall system
c
         if (debug) then
            write (iout,20)  (vang(i),i=1,3),erot
   20       format (' System Angular Velocity : ',3d12.2,
     &              /,' Rotational Kinetic Energy :',13x,f12.4,
     &                 ' Kcal/mole')
         end if
      end if
c
c     perform deallocation of some local arrays
c
      if (.not.use_bounds .and. integrate.eq.'RIGIDBODY') then
         deallocate (xcm)
         deallocate (ycm)
         deallocate (zcm)
      end if
      return
      end
