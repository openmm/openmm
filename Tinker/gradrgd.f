c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine gradrgd  --  energy & gradient of rigid body  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "gradrgd" calls subroutines to calculate the potential energy
c     and first derivatives with respect to rigid body coordinates
c
c
      subroutine gradrgd (energy,derivs)
      use sizes
      use atoms
      use group
      use rigid
      implicit none
      integer i,j,k
      integer init,stop
      real*8 energy
      real*8 xcm,ycm,zcm
      real*8 xterm,yterm,zterm
      real*8 phi,cphi,sphi
      real*8 theta,ctheta,stheta
      real*8 ephi(3),etheta(3)
      real*8 epsi(3),tau(3)
      real*8 derivs(6,*)
      real*8, allocatable :: g(:,:)
c
c
c     zero out the total of rigid body derivative components
c
      do i = 1, ngrp
         do j = 1, 6
            derivs(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (g(3,n))
c
c     calculate the energy and Cartesian first derivatives
c
      call gradient (energy,g)
c
c     compute the rigid body gradient components for each group
c
      do i = 1, ngrp
         init = igrp(1,i)
         stop = igrp(2,i)
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
c
c     get unit vectors along the phi, theta and psi rotation axes
c
         ephi(1) = 0.0d0
         ephi(2) = 0.0d0
         ephi(3) = 1.0d0
         etheta(1) = -sphi
         etheta(2) = cphi
         etheta(3) = 0.0d0
         epsi(1) = ctheta * cphi
         epsi(2) = ctheta * sphi
         epsi(3) = -stheta
c
c     find the rigid body gradients for translations
c
         do j = init, stop
            k = kgrp(j)
            derivs(1,i) = derivs(1,i) + g(1,k)
            derivs(2,i) = derivs(2,i) + g(2,k)
            derivs(3,i) = derivs(3,i) + g(3,k)
         end do
c
c     accumulate the moment arm along each axis of rotation
c
         do j = 1, 3
            tau(j) = 0.0d0
         end do
         do j = init, stop
            k = kgrp(j)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            tau(1) = tau(1) + yterm*g(3,k) - zterm*g(2,k)
            tau(2) = tau(2) + zterm*g(1,k) - xterm*g(3,k)
            tau(3) = tau(3) + xterm*g(2,k) - yterm*g(1,k)
         end do
c
c     find the rigid body gradients for rotations
c
         do j = 1, 3
            derivs(4,i) = derivs(4,i) + tau(j)*ephi(j)
            derivs(5,i) = derivs(5,i) + tau(j)*etheta(j)
            derivs(6,i) = derivs(6,i) + tau(j)*epsi(j)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (g)
      return
      end
