c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine impose  --  superimpose two coordinate sets  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "impose" performs the least squares best superposition
c     of two atomic coordinate sets via a quaternion method;
c     upon return, the first coordinate set is unchanged while
c     the second set is translated and rotated to give best fit;
c     the final root mean square fit is returned in "rmsvalue"
c
c
      subroutine impose (n1,x1,y1,z1,n2,x2,y2,z2,rmsvalue)
      use sizes
      use align
      use inform
      use iounit
      implicit none
      integer i,n1,n2,nmax
      real*8 xmid,ymid,zmid
      real*8 rmsvalue,rmsfit
      real*8 x1(*),x2(*)
      real*8 y1(*),y2(*)
      real*8 z1(*),z2(*)
c
c
c     perform dynamic allocation of some global arrays
c
      nmax = max(n1,n2)
      if (.not. allocated(ifit))  allocate (ifit(2,nmax))
      if (.not. allocated(wfit))  allocate (wfit(nmax))
c
c     superimpose the full structures if not specified
c
      if (nfit .eq. 0) then
         nfit = min(n1,n2)
         do i = 1, nfit
            ifit(1,i) = i
            ifit(2,i) = i
            wfit(i) = 1.0d0
         end do
      end if
c
c     if the weights are all zero, set them to unity
c
      do i = 1, nfit
         if (wfit(i) .ne. 0.0d0)  goto 10
      end do
      do i = 1, nfit
         wfit(i) = 1.0d0
      end do
   10 continue
c
c     find the rms fit of input coordinates
c
      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
         write (iout,20)  rmsvalue
   20    format (/,' IMPOSE  --  Input Coordinates',12x,f12.6)
      end if
c
c     superimpose the centroids of active atom pairs
c
      call center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid)
      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
         write (iout,30)  rmsvalue
   30    format (' IMPOSE  --  After Translation',12x,f12.6)
      end if
c
c     use a quaternion method to achieve the superposition
c
      call quatfit (n1,x1,y1,z1,n2,x2,y2,z2)
      rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
      if (verbose) then
         write (iout,40)  rmsvalue
   40    format (' IMPOSE  --  After Rotation',15x,f12.6)
      end if
c
c     translate both coordinate sets so as to return
c     the first set to its original position
c
      do i = 1, n1
         x1(i) = x1(i) + xmid
         y1(i) = y1(i) + ymid
         z1(i) = z1(i) + zmid
      end do
      do i = 1, n2
         x2(i) = x2(i) + xmid
         y2(i) = y2(i) + ymid
         z2(i) = z2(i) + zmid
      end do
      return
      end
