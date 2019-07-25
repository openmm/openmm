c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  function erf  --  evaluate the standard error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "erf" computes a numerical approximation to the value of
c     the error function via a Chebyshev approximation
c
c
      function erf (x)
      implicit none
      integer mode
      real*8 erf,x
      real*8 result
c
c
c     compute the error function via Chebyshev fitting
c
      mode = 0
      call erfcore (x,result,mode)
      erf = result
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function erfc  --  evaluate complementary error function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erfc" computes a numerical approximation to the value of the
c     complementary error function via a Chebyshev approximation
c
c
      function erfc (x)
      implicit none
      integer mode
      real*8 erfc,x
      real*8 result
c
c
c     get the complementary error function via Chebyshev fitting
c
      mode = 1
      call erfcore (x,result,mode)
      erfc = result
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erfcore  --  erf and erfc via Chebyshev approx  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erfcore" evaluates erf(x) or erfc(x) for a real argument x;
c     when called with mode set to 0 it returns erf, a mode of 1
c     returns erfc; uses rational functions that approximate erf(x)
c     and erfc(x) to at least 18 significant decimal digits
c
c     literature reference:
c
c     W. J. Cody, "Rational Chebyshev Approximations for the Error
c     Function", Mathematics of Computation, 631-638, 1969
c
c     adapted from an original program written by W. J. Cody,
c     Mathematics and Computer Science Division, Argonne National
c     Laboratory, Argonne, IL 60439
c
c     machine-dependent constants:
c
c     xtiny   argument below which erf(x) may be represented by
c             2*x/sqrt(pi) and above which x*x won't underflow;
c             a conservative value is the largest machine number
c             X such that 1.0 + X = 1.0 to machine precision
c
c     xbig    largest argument acceptable for erfc; solution to
c             the equation:  W(x) * (1-0.5/x**2) = XMIN, where
c             W(x) = exp(-x*x)/[x*sqrt(pi)]
c
c
      subroutine erfcore (arg,result,mode)
      implicit none
      integer i,mode
      real*8 arg,result
      real*8 x,y,ysq
      real*8 del,sqrpi
      real*8 xnum,xden
      real*8 xtiny,xbig
      real*8 a(5),b(4)
      real*8 c(9),d(8)
      real*8 p(6),q(5)
c
c     mathematical and machine-dependent constants
c
      data sqrpi  / 5.6418958354775628695d-1 /
      data xtiny  / 1.11d-16 /
      data xbig   / 26.543d0 /
c
c     coefficients for approximation to erf in first interval
c
      data a  / 3.16112374387056560d0,  1.13864154151050156d2,
     &          3.77485237685302021d2,  3.20937758913846947d3,
     &          1.85777706184603153d-1 /
      data b  / 2.36012909523441209d1,  2.44024637934444173d2,
     &          1.28261652607737228d3,  2.84423683343917062d3 /
c
c     coefficients for approximation to erfc in second interval
c
      data c  / 5.64188496988670089d-1, 8.88314979438837594d0,
     &          6.61191906371416295d1,  2.98635138197400131d2,
     &          8.81952221241769090d2,  1.71204761263407058d3,
     &          2.05107837782607147d3,  1.23033935479799725d3,
     &          2.15311535474403846d-8 /
      data d  / 1.57449261107098347d1,  1.17693950891312499d2,
     &          5.37181101862009858d2,  1.62138957456669019d3,
     &          3.29079923573345963d3,  4.36261909014324716d3,
     &          3.43936767414372164d3,  1.23033935480374942d3 /
c
c     coefficients for approximation to erfc in third interval
c
      data p  / 3.05326634961232344d-1, 3.60344899949804439d-1,
     &          1.25781726111229246d-1, 1.60837851487422766d-2,
     &          6.58749161529837803d-4, 1.63153871373020978d-2 /
      data q  / 2.56852019228982242d0,  1.87295284992346047d0,
     &          5.27905102951428412d-1, 6.05183413124413191d-2,
     &          2.33520497626869185d-3 /
c
c
c     store the argument and its absolute value
c
      x = arg
      y = abs(x)
c
c     evaluate error function for |x| less than 0.46875
c
      if (y .le. 0.46875d0) then
         ysq = 0.0d0
         if (y .gt. xtiny)  ysq = y * y
         xnum = a(5) * ysq
         xden = ysq
         do i = 1, 3
            xnum = (xnum + a(i)) * ysq
            xden = (xden + b(i)) * ysq
         end do
         result = x * (xnum + a(4)) / (xden + b(4))
         if (mode .ne. 0)  result = 1.0d0 - result
c
c     get complementary error function for 0.46875 <= |x| <= 4.0
c
      else if (y .le. 4.0d0) then
         xnum = c(9) * y
         xden = y
         do i = 1, 7
            xnum = (xnum + c(i)) * y
            xden = (xden + d(i)) * y
         end do
         result = (xnum + c(8)) / (xden + d(8))
         ysq = aint(16.0d0*y) / 16.0d0
         del = (y-ysq) * (y+ysq)
c        result = exp(-ysq*ysq) * exp(-del) * result
         result = exp(-ysq*ysq-del) * result
         if (mode .eq. 0) then
            result = 1.0d0 - result
            if (x .lt. 0.0d0)  result = -result
         else
            if (x .lt. 0.0d0)  result = 2.0d0 - result
         end if
c
c     get complementary error function for |x| greater than 4.0
c
      else
         result = 0.0d0
         if (y .lt. xbig) then
            ysq = 1.0d0 / (y * y)
            xnum = p(6) * ysq
            xden = ysq
            do i = 1, 4
               xnum = (xnum + p(i)) * ysq
               xden = (xden + q(i)) * ysq
            end do
            result = ysq * (xnum + p(5)) / (xden + q(5))
            result = (sqrpi - result) / y
            ysq = aint(16.0d0*y) / 16.0d0
            del = (y-ysq) * (y+ysq)
c           result = exp(-ysq*ysq) * exp(-del) * result
            result = exp(-ysq*ysq-del) * result
         end if
         if (mode .eq. 0) then
            result = 1.0d0 - result
            if (x .lt. 0.0d0)  result = -result
         else
            if (x .lt. 0.0d0)  result = 2.0d0 - result
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function erfinv  --  evaluate the error function inverse  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erfinv" evaluates the inverse of the error function for
c     an argument in the range (-1,1) using a rational function
c     approximation followed by cycles of Newton-Raphson correction
c
c     adapted from the pseudocode for the Matlab function of the
c     same name; Matlab, version 4.2c, March 1995
c
c
      function erfinv (x)
      use iounit
      use math
      implicit none
      real*8 erfinv,erf
      real*8 x,y,z
      real*8 a(4),b(4)
      real*8 c(4),d(2)
      external erf
c
c     coefficients for approximation to erfinv in central range
c
      data a  /  0.886226899d0, -1.645349621d0,
     &           0.914624893d0, -0.140543331d0 /
      data b  / -2.118377725d0,  1.442710462d0,
     &          -0.329097515d0,  0.012229801d0 /
c
c     coefficients for approximation to erfinv near endpoints
c
      data c  / -1.970840454d0, -1.624906493d0,
     &           3.429567803d0,  1.641345311d0 /
      data d  /  3.543889200d0,  1.637067800d0 /
c
c
c     get an initial estimate for the inverse error function
c
      if (abs(x) .le. 0.7d0) then
         y = x * x
         z = x * (((a(4)*y+a(3))*y+a(2))*y+a(1))
     &              / ((((b(4)*y+b(3))*y+b(2))*y+b(1))*y+1.0d0)
      else if (x.gt.0.7d0 .and. x.lt.1.0d0) then
         y = sqrt(-log((1.0d0-x)/2.0d0))
         z = (((c(4)*y+c(3))*y+c(2))*y+c(1)) / ((d(2)*y+d(1))*y+1.0d0)
      else if (x.lt.-0.7d0 .and. x.gt.-1.0d0) then
         y = sqrt(-log((1.0d0+x)/2.0d0))
         z = -(((c(4)*y+c(3))*y+c(2))*y+c(1)) / ((d(2)*y+d(1))*y+1.0d0)
      else
         write (iout,10)
   10    format (/,' ERFINV  --  Illegal Argument to Inverse',
     &              ' Error Function')
         call fatal
      end if
c
c     use two steps of Newton-Raphson correction to increase accuracy
c
      z = z - (erf(z) - x) / (2.0d0/sqrtpi * exp(-z*z))
      z = z - (erf(z) - x) / (2.0d0/sqrtpi * exp(-z*z))
      erfinv = z
      return
      end
