c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine cspline  --  periodic interpolating cube spline  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "cspline" computes the coefficients for a periodic interpolating
c     cubic spline
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 10.1.2  [see routine "isplpe"]
c
c
      subroutine cspline (n,xn,fn,b,c,d,h,du,dm,rc,rs)
      use iounit
      implicit none
      integer i,n,iflag
      real*8 eps,average
      real*8 temp1,temp2
      real*8 xn(0:*)
      real*8 fn(0:*)
      real*8 b(0:*)
      real*8 c(0:*)
      real*8 d(0:*)
      real*8 h(0:*)
      real*8 du(0:*)
      real*8 dm(0:*)
      real*8 rc(0:*)
      real*8 rs(0:*)
c
c
c     check the periodicity of fn, and for subsequent call
c
      eps = 0.000001d0
      if (abs(fn(n)-fn(0)) .gt. eps) then
         write (iout,10)  fn(0),fn(n)
   10    format (/,' CSPLINE  --  Warning, Non-Periodic Input',
     &              ' Values',2f12.5)
      end if
      average = 0.5d0 * (fn(0) + fn(n))
      fn(0) = average
      fn(n) = average
c
c     get auxiliary variables and matrix elements on first call
c
      do i = 0, n-1
         h(i) = xn(i+1) - xn(i)
      end do
      h(n) = h(0)
      do i = 1, n-1
         du(i) = h(i)
      end do
      du(n) = h(0)
      do i = 1, n
         dm(i) = 2.0d0 * (h(i-1)+h(i))
      end do
c
c     compute the right hand side
c
      temp1 = (fn(1)-fn(0)) / h(0)
      do i = 1, n-1, 1
         temp2 = (fn(i+1)-fn(i)) / h(i)
         rs(i)  = 3.0d0 * (temp2-temp1)
         temp1 = temp2
      end do
      rs(n) = 3.0d0 * ((fn(1)-fn(0))/h(0)-temp1)
c
c     solve the linear system with factorization
c
      call cytsy (n,dm,du,rc,rs,c,iflag)
      if (iflag .ne. 1)  return
c
c     compute remaining spline coefficients
c
      c(0) = c(n)
      do i = 0, n-1
         b(i) = (fn(i+1)-fn(i))/h(i) - h(i)/3.0d0*(c(i+1)+2.0d0*c(i))
         d(i) = (c(i+1)-c(i)) / (3.0d0*h(i))
      end do
      b(n) = (fn(1)-fn(n))/h(n) - h(n)/3.0d0*(c(1)+2.0d0*c(n))
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine cytsy  --  solve cyclic tridiagonal system  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "cytsy" solves a system of linear equations for a cyclically
c     tridiagonal, symmetric, positive definite matrix
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsy (n,dm,du,cr,rs,x,iflag)
      implicit none
      integer n,iflag
      real*8 dm(0:*)
      real*8 du(0:*)
      real*8 cr(0:*)
      real*8 rs(0:*)
      real*8 x(0:*)
c
c
c     factorization of the input matrix
c
      iflag = -2
      if (n .lt. 3)  return
      call cytsyp (n,dm,du,cr,iflag)
c
c     update and back substitute as necessary
c
      if (iflag .eq. 1)  call cytsys (n,dm,du,cr,rs,x)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cytsyp  --  tridiagonal Cholesky factorization  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cytsyp" finds the Cholesky factors of a cyclically tridiagonal
c     symmetric, positive definite matrix given by two vectors
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsyp (n,dm,du,cr,iflag)
      implicit none
      integer i,n,iflag
      real*8 eps,row,d
      real*8 temp1,temp2
      real*8 dm(0:*)
      real*8 du(0:*)
      real*8 cr(0:*)
c
c
c     set error bound and test for condition n greater than 2
c
      eps = 0.00000001d0
      iflag = -2
      if (n .lt. 3)  return
c
c     checking to see if matrix is positive definite
c
      row = abs(dm(1)) + abs(du(1)) + abs(du(n))
      if (row .eq. 0.0d0) then
         iflag = 0
         return
      end if
      d = 1.0d0 / row
      if (dm(1) .lt. 0.0d0) then
         iflag = -1
         return
      else if (abs(dm(1))*d .le. eps) then
         iflag = 0
         return
      end if
c
c     factoring a while checking for a positive definite and strong
c     nonsingular matrix a
c
      temp1 = du(1)
      du(1) = du(1) / dm(1)
      cr(1) = du(n) / dm(1)
      do i = 2, n-1
         row = abs(dm(i)) + abs(du(i)) + abs(temp1)
         if (row .eq. 0.0d0) then
            iflag = 0
            return
         end if
         d = 1.0d0 / row
         dm(i) = dm(i) - temp1*du(i-1)
         if (dm(i) .lt. 0.0d0) then
            iflag = -1
            return
         else if (abs(dm(i))*d .le. eps) then
            iflag = 0
            return
         end if
         if (i .lt. (n-1)) then
            cr(i) = -temp1 * cr(i-1) / dm(i)
            temp1 = du(i)
            du(i) = du(i) / dm(i)
         else
            temp2 = du(i)
            du(i) = (du(i) - temp1*cr(i-1)) / dm(i)
         end if
      end do
      row = abs(du(n)) + abs(dm(n)) + abs(temp2)
      if (row .eq. 0.0d0) then
         iflag = 0
         return
      end if
      d = 1.0d0 / row
      dm(n) = dm(n) - dm(n-1)*du(n-1)*du(n-1)
      temp1 = 0.0d0
      do i = 1, n-2
         temp1 = temp1 + dm(i)*cr(i)*cr(i)
      end do
      dm(n) = dm(n) - temp1
      if (dm(n) .lt. 0) then
         iflag = -1
         return
      else if (abs(dm(n))*d .le. eps) then
         iflag = 0
         return
      end if
      iflag = 1
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cytsys  --  tridiagonal solution from factors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cytsys" solves a cyclically tridiagonal linear system
c     given the Cholesky factors
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsys (n,dm,du,cr,rs,x)
      implicit none
      integer i,n
      real*8 sum,temp
      real*8 dm(0:*)
      real*8 du(0:*)
      real*8 cr(0:*)
      real*8 rs(0:*)
      real*8 x(0:*)
c
c
c     updating phase
c
      temp = rs(1)
      rs(1) = temp / dm(1)
      sum = cr(1) * temp
      do i = 2, n-1
         temp = rs(i) - du(i-1)*temp
         rs(i) = temp / dm(i)
         if (i .ne. (n-1))  sum = sum + cr(i)*temp
      end do
      temp = rs(n) - du(n-1)*temp
      temp = temp - sum
      rs(n) = temp / dm(n)
c
c     back substitution phase
c
      x(n) = rs(n)
      x(n-1) = rs(n-1) - du(n-1)*x(n)
      do i = n-2, 1, -1
         x(i) = rs(i) - du(i)*x(i+1) - cr(i)*x(n)
      end do
      return
      end
