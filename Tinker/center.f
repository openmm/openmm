c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine center  --  superimpose structure centroids  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "center" moves the weighted centroid of each coordinate
c     set to the origin during least squares superposition
c
c
      subroutine center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid)
      use align
      implicit none
      integer i,k,n1,n2
      real*8 weigh,norm
      real*8 xmid,ymid,zmid
      real*8 x1(*),x2(*)
      real*8 y1(*),y2(*)
      real*8 z1(*),z2(*)
c
c
c     find the weighted centroid of the second
c     structure and translate it to the origin
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(2,i)
         weigh = wfit(i)
         xmid = xmid + x2(k)*weigh
         ymid = ymid + y2(k)*weigh
         zmid = zmid + z2(k)*weigh
         norm = norm + weigh
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n2
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do
c
c     now repeat for the first structure, note
c     that this centroid position gets returned
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(1,i)
         weigh = wfit(i)
         xmid = xmid + x1(k)*weigh
         ymid = ymid + y1(k)*weigh
         zmid = zmid + z1(k)*weigh
         norm = norm + weigh
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n1
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do
      return
      end
