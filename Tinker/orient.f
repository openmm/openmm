c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine orient  --  rigid body reference coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "orient" computes a set of reference Cartesian coordinates
c     in standard orientation for each rigid body atom group
c
c
      subroutine orient
      use sizes
      use atoms
      use group
      use rigid
      implicit none
      integer i,j,k
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 xterm,yterm,zterm
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3)
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(xrb))  allocate (xrb(n))
      if (.not. allocated(yrb))  allocate (yrb(n))
      if (.not. allocated(zrb))  allocate (zrb(n))
      if (.not. allocated(rbc))  allocate (rbc(6,ngrp))
c
c     use current coordinates as default reference coordinates
c
      do i = 1, n
         xrb(i) = x(i)
         yrb(i) = y(i)
         zrb(i) = z(i)
      end do
c
c     compute the rigid body coordinates for each atom group
c
      call xyzrigid
c
c     get the center of mass and Euler angles for each group
c
      do i = 1, ngrp
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         psi = rbc(6,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
c
c     construct the rotation matrix from Euler angle values
c
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
c
c     translate and rotate each atom group into inertial frame
c
         init = igrp(1,i)
         stop = igrp(2,i)
         do j = init, stop
            k = kgrp(j)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            xrb(k) = a(1,1)*xterm + a(1,2)*yterm + a(1,3)*zterm
            yrb(k) = a(2,1)*xterm + a(2,2)*yterm + a(2,3)*zterm
            zrb(k) = a(3,1)*xterm + a(3,2)*yterm + a(3,3)*zterm
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine xyzrigid  --  determine rigid body coordinates  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xyzrigid" computes the center of mass and Euler angle rigid
c     body coordinates for each atom group in the system
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine xyzrigid
      use sizes
      use atoms
      use atomid
      use group
      use rigid
      implicit none
      integer i,j,k,m
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 weigh,total,dot
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 moment(3)
      real*8 vec(3,3)
      real*8 tensor(3,3)
      real*8 a(3,3)
c
c
c     get the first and last atom in the current group
c
      do i = 1, ngrp
         init = igrp(1,i)
         stop = igrp(2,i)
c
c     compute the position of the group center of mass
c
         total = 0.0d0
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do j = init, stop
            k = kgrp(j)
            weigh = mass(k)
            total = total + weigh
            xcm = xcm + x(k)*weigh
            ycm = ycm + y(k)*weigh
            zcm = zcm + z(k)*weigh
         end do
         if (total .ne. 0.0d0) then
            xcm = xcm / total
            ycm = ycm / total
            zcm = zcm / total
         end if
c
c     compute and then diagonalize the inertia tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         do j = init, stop
            k = kgrp(j)
            weigh = mass(k)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            xx = xx + xterm*xterm*weigh
            xy = xy + xterm*yterm*weigh
            xz = xz + xterm*zterm*weigh
            yy = yy + yterm*yterm*weigh
            yz = yz + yterm*zterm*weigh
            zz = zz + zterm*zterm*weigh
         end do
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
         call jacobi (3,tensor,moment,vec)
c
c     select the direction for each principle moment axis
c
         do m = 1, 2
            do j = init, stop
               k = kgrp(j)
               xterm = vec(1,m) * (x(k)-xcm)
               yterm = vec(2,m) * (y(k)-ycm)
               zterm = vec(3,m) * (z(k)-zcm)
               dot = xterm + yterm + zterm
               if (dot .lt. 0.0d0) then
                  vec(1,m) = -vec(1,m)
                  vec(2,m) = -vec(2,m)
                  vec(3,m) = -vec(3,m)
               end if
               if (dot .ne. 0.0d0)  goto 10
            end do
   10       continue
         end do
c
c     moment axes must give a right-handed coordinate system
c
         xterm = vec(1,1) * (vec(2,2)*vec(3,3)-vec(2,3)*vec(3,2))
         yterm = vec(2,1) * (vec(1,3)*vec(3,2)-vec(1,2)*vec(3,3))
         zterm = vec(3,1) * (vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2))
         dot = xterm + yterm + zterm
         if (dot .lt. 0.0d0) then
            do j = 1, 3
               vec(j,3) = -vec(j,3)
            end do
         end if
c
c     principal moment axes form rows of Euler rotation matrix
c
         do k = 1, 3
            do j = 1, 3
               a(k,j) = vec(j,k)
            end do
         end do
c
c     compute Euler angles consistent with the rotation matrix
c
         call roteuler (a,phi,theta,psi)
c
c     set the rigid body coordinates for each atom group
c
         rbc(1,i) = xcm
         rbc(2,i) = ycm
         rbc(3,i) = zcm
         rbc(4,i) = phi
         rbc(5,i) = theta
         rbc(6,i) = psi
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine roteuler  --  rotation matrix to Euler angles   ##
c     ##                                                             ##
c     #################################################################
c
c
c     "roteuler" computes a set of Euler angle values consistent
c     with an input rotation matrix
c
c
      subroutine roteuler (a,phi,theta,psi)
      use math
      implicit none
      integer i
      real*8 phi,theta,psi,eps
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3),b(3)
      logical flip(3)
c
c
c     set the tolerance for Euler angles and rotation elements
c
      eps = 1.0d-7
c
c     get a trial value of theta from a single rotation element
c
      theta = asin(min(1.0d0,max(-1.0d0,-a(1,3))))
      ctheta = cos(theta)
      stheta = -a(1,3)
c
c     set the phi/psi difference when theta is either 90 or -90
c
      if (abs(ctheta) .le. eps) then
         phi = 0.0d0
         if (abs(a(3,1)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,-a(2,1)/a(1,3))))
         else if (abs(a(2,1)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,-a(3,1)/a(1,3))))
         else
            psi = atan(a(2,1)/a(3,1))
         end if
c
c     set the phi and psi values for all other theta values
c
      else
         if (abs(a(1,1)) .lt. eps) then
            phi = asin(min(1.0d0,max(-1.0d0,a(1,2)/ctheta)))
         else if (abs(a(1,2)) .lt. eps) then
            phi = acos(min(1.0d0,max(-1.0d0,a(1,1)/ctheta)))
         else
            phi = atan(a(1,2)/a(1,1))
         end if
         if (abs(a(3,3)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,a(2,3)/ctheta)))
         else if (abs(a(2,3)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,a(3,3)/ctheta)))
         else
            psi = atan(a(2,3)/a(3,3))
         end if
      end if
c
c     find sine and cosine of the trial phi and psi values
c
      cphi = cos(phi)
      sphi = sin(phi)
      cpsi = cos(psi)
      spsi = sin(psi)
c
c     reconstruct the diagonal of the rotation matrix
c
      b(1) = ctheta * cphi
      b(2) = spsi*stheta*sphi + cpsi*cphi
      b(3) = ctheta * cpsi
c
c     compare the correct matrix diagonal to rebuilt diagonal
c
      do i = 1, 3
         flip(i) = .false.
         if (abs(a(i,i)-b(i)) .gt. eps)  flip(i) = .true.
      end do
c
c     alter Euler angles to get correct rotation matrix values
c
      if (flip(1) .and. flip(2))  phi = phi - sign(pi,phi)
      if (flip(1) .and. flip(3))  theta = -theta + sign(pi,theta)
      if (flip(2) .and. flip(3))  psi = psi - sign(pi,psi)
c
c     convert maximum negative angles to positive values
c
      if (phi .le. -pi)  phi = pi
      if (theta .le. -pi)  theta = pi
      if (psi .le. -pi)  psi = pi
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rigidxyz  --  rigid body to Cartesian coords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rigidxyz" computes Cartesian coordinates for a rigid body
c     group via rotation and translation of reference coordinates
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine rigidxyz
      use sizes
      use atoms
      use group
      use rigid
      implicit none
      integer i,j,k
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 xterm,yterm,zterm
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3)
c
c
c     get the center of mass and Euler angles for each group
c
      do i = 1, ngrp
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         psi = rbc(6,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
c
c     construct the rotation matrix from Euler angle values
c
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
c
c     rotate and translate reference coordinates into global frame
c
         init = igrp(1,i)
         stop = igrp(2,i)
         do j = init, stop
            k = kgrp(j)
            xterm = xrb(k)
            yterm = yrb(k)
            zterm = zrb(k)
            x(k) = a(1,1)*xterm + a(2,1)*yterm + a(3,1)*zterm + xcm
            y(k) = a(1,2)*xterm + a(2,2)*yterm + a(3,2)*zterm + ycm
            z(k) = a(1,3)*xterm + a(2,3)*yterm + a(3,3)*zterm + zcm
         end do
      end do
      return
      end
