c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine gyrate  --  compute the radius of gyration  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "gyrate" computes the radius of gyration of a molecular system
c     from its atomic coordinates; only active atoms are included
c
c
      subroutine gyrate (rg)
      use sizes
      use atoms
      use usage
      implicit none
      integer i,k
      real*8 rg,xc,yc,zc
c
c
c     find the centroid of the atomic coordinates
c
      xc = 0.0d0
      yc = 0.0d0
      zc = 0.0d0
      do i = 1, nuse
         k = iuse(i)
         xc = xc + x(k)
         yc = yc + y(k)
         zc = zc + z(k)
      end do
      xc = xc / dble(nuse)
      yc = yc / dble(nuse)
      zc = zc / dble(nuse)
c
c     compute and print out the radius of gyration
c
      rg = 0.0d0
      do i = 1, nuse
         k = iuse(i)
         rg = rg + (x(k)-xc)**2 + (y(k)-yc)**2 + (z(k)-zc)**2
      end do
      rg = sqrt(rg/dble(nuse))
      return
      end
