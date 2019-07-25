c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  function rmsfit  --  rms deviation for paired atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "rmsfit" computes the rms fit of two coordinate sets
c
c
      function rmsfit (x1,y1,z1,x2,y2,z2)
      use sizes
      use align
      implicit none
      integer i,i1,i2
      real*8 rmsfit,rmsterm
      real*8 xr,yr,zr,dist2
      real*8 weigh,norm
      real*8 x1(*),x2(*)
      real*8 y1(*),y2(*)
      real*8 z1(*),z2(*)
c
c
c     compute the rms fit over superimposed atom pairs
c
      rmsfit = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weigh = wfit(i)
         xr = x1(i1) - x2(i2)
         yr = y1(i1) - y2(i2)
         zr = z1(i1) - z2(i2)
         dist2 = xr**2 + yr**2 + zr**2
         norm = norm + weigh
         rmsterm = dist2 * weigh
         rmsfit = rmsfit + rmsterm
      end do
      rmsfit = sqrt(rmsfit/norm)
      return
      end
