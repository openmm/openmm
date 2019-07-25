c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Chuanjie Wu and Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##   subroutine nspline  --  nonperiodic natural cubic spline   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "nspline" computes coefficients for an nonperiodic cubic spline
c     with natural boundary conditions where the first and last second
c     derivatives are already known
c
c
      subroutine nspline (n,x0,y0,s1,s2,h,g,dy,dla,dmu)
      implicit none
      integer i,n
      real*8 t,y21,y2n
      real*8 x0(0:*)
      real*8 y0(0:*)
      real*8 s1(0:*)
      real*8 s2(0:*)
      real*8 h(0:*)
      real*8 g(0:*)
      real*8 dy(0:*)
      real*8 dla(0:*)
      real*8 dmu(0:*)
c
c
c     set first and last second deriviatives to zero
c
      y21 = 0.0d0
      y2n = 0.0d0
c
c     find the intervals to be used
c
      do i = 0, n-1
         h(i) = x0(i+1) - x0(i)
         dy(i) = (y0(i+1)-y0(i)) / h(i)
      end do
c
c     calculate the spline coeffcients
c
      do i = 1, n-1
         dla(i) = h(i) / (h(i)+h(i-1))
         dmu(i) = 1.0d0 - dla(i)
         g(i) = 3.0d0 * (dla(i)*dy(i-1)+dmu(i)*dy(i))
      end do
c
c     set the initial value via natural boundary condition
c
      dla(n) = 1.0d0
      dla(0) = 0.0d0
      dmu(n) = 0.0d0
      dmu(0) = 1.0d0
      g(0) = 3.0d0*dy(0) - 0.5d0*h(0)*y21
      g(n) = 3.0d0*dy(n-1) + 0.5d0*h(n-1)*y2n
c
c     solve the triagonal system of linear equations
c
      dmu(0) = 0.5d0 * dmu(0)
      g(0) = 0.5d0 * g(0)
      do i = 1, n
         t = 2.0d0 - dmu(i-1)*dla(i)
         dmu(i) = dmu(i) / t
         g(i) = (g(i)-g(i-1)*dla(i)) / t
      end do
      do i = n-1, 0, -1
         g(i) = g(i) - dmu(i)*g(i+1)
      end do
c
c     get the first derivative at each grid point
c
      do i = 0, n
         s1(i) = g(i)
      end do
c
c     get the second derivative at each grid point
c
      s2(0) = y21
      s2(n) = y2n
      do i = 1, n-1
         s2(i) = 6.0d0*(y0(i+1)-y0(i))/(h(i)*h(i))
     &              - 4.0d0*s1(i)/h(i) - 2.0d0*s1(i+1)/h(i)
      end do
      return
      end
