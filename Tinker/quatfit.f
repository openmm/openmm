c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine quatfit  --  quaternion superposition of coords  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "quatfit" uses a quaternion-based method to achieve the best
c     fit superposition of two sets of coordinates
c
c     literature reference:
c
c     S. K. Kearsley, "On the Orthogonal Transformation Used for
c     Structural Comparisons", Acta Crystallographica Section A,
c     45, 208-210 (1989)
c
c     adapted from an original program written by D. J. Heisterberg,
c     Ohio Supercomputer Center, Columbus, OH
c
c
      subroutine quatfit (n1,x1,y1,z1,n2,x2,y2,z2)
      use sizes
      use align
      implicit none
      integer i,i1,i2,n1,n2
      real*8 weigh,xrot,yrot,zrot
      real*8 xxyx,xxyy,xxyz
      real*8 xyyx,xyyy,xyyz
      real*8 xzyx,xzyy,xzyz
      real*8 q(4),d(4)
      real*8 x1(*),x2(*)
      real*8 y1(*),y2(*)
      real*8 z1(*),z2(*)
      real*8 rot(3,3)
      real*8 c(4,4),v(4,4)
c
c
c     build the upper triangle of the quadratic form matrix
c
      xxyx = 0.0d0
      xxyy = 0.0d0
      xxyz = 0.0d0
      xyyx = 0.0d0
      xyyy = 0.0d0
      xyyz = 0.0d0
      xzyx = 0.0d0
      xzyy = 0.0d0
      xzyz = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weigh = wfit(i)
         xxyx = xxyx + weigh*x1(i1)*x2(i2)
         xxyy = xxyy + weigh*y1(i1)*x2(i2)
         xxyz = xxyz + weigh*z1(i1)*x2(i2)
         xyyx = xyyx + weigh*x1(i1)*y2(i2)
         xyyy = xyyy + weigh*y1(i1)*y2(i2)
         xyyz = xyyz + weigh*z1(i1)*y2(i2)
         xzyx = xzyx + weigh*x1(i1)*z2(i2)
         xzyy = xzyy + weigh*y1(i1)*z2(i2)
         xzyz = xzyz + weigh*z1(i1)*z2(i2)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy
c
c     diagonalize the quadratic form matrix
c
      call jacobi (4,c,d,v)
c
c     extract the desired quaternion
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble rotation matrix that superimposes the molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
      do i = 1, n2
         xrot = x2(i)*rot(1,1) + y2(i)*rot(1,2) + z2(i)*rot(1,3)
         yrot = x2(i)*rot(2,1) + y2(i)*rot(2,2) + z2(i)*rot(2,3)
         zrot = x2(i)*rot(3,1) + y2(i)*rot(3,2) + z2(i)*rot(3,3)
         x2(i) = xrot
         y2(i) = yrot
         z2(i) = zrot
      end do
      return
      end
