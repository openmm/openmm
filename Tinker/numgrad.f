c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine numgrad  --  numerical gradient of a function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "numgrad" computes the gradient of the objective function
c     "fvalue" with respect to Cartesian coordinates of the atoms
c     via a one-sided or two-sided numerical differentiation
c
c
      subroutine numgrad (fvalue,g,eps)
      use sizes
      use atoms
      implicit none
      integer i
      real*8 fvalue,eps
      real*8 f,f0,old
      real*8 g(3,*)
      logical twosided
      external fvalue
c
c
c     chose between use of one-sided or two-sided gradient
c
      twosided = .true.
      if (.not. twosided)  f0 = fvalue ()
c
c     compute the numerical gradient from function values
c
      do i = 1, n
         old = x(i)
         if (twosided) then
            x(i) = x(i) - 0.5d0*eps
            f0 = fvalue ()
         end if
         x(i) = x(i) + eps
         f = fvalue ()
         x(i) = old
         g(1,i) = (f - f0) / eps
         old = y(i)
         if (twosided) then
            y(i) = y(i) - 0.5d0*eps
            f0 = fvalue ()
         end if
         y(i) = y(i) + eps
         f = fvalue ()
         y(i) = old
         g(2,i) = (f - f0) / eps
         old = z(i)
         if (twosided) then
            z(i) = z(i) - 0.5d0*eps
            f0 = fvalue ()
         end if
         z(i) = z(i) + eps
         f = fvalue ()
         z(i) = old
         g(3,i) = (f - f0) / eps
      end do
      return
      end
