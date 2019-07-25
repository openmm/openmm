c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine settime  --  initialize wall clock and CPU time  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "settime" initializes the wall clock and elapsed CPU times
c
c
      subroutine settime
      use chrono
      implicit none
      integer count,crate,cmax
      real time
c
c
c     set wall clock and cpu time via Fortran intrinsic functions
c
      call system_clock (count,crate,cmax)
      twall = dble(count) / dble(crate)
      call cpu_time (time)
      tcpu = dble(time)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gettime  --  elapsed wall clock and CPU times  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gettime" finds the elapsed wall clock and CPU times in seconds
c     since the last call to "settime"
c
c
      subroutine gettime (wall,cpu)
      use chrono
      implicit none
      integer count,crate,cmax
      real time
      real*8 wall,cpu
      real*8 twall0,tcpu0
c
c
c     get total wall clock time via Fortran intrinsic functions
c
      twall0 = twall
      call system_clock (count,crate,cmax)
      twall = dble(count) / dble(crate)
      wall = twall - twall0
      if (wall .lt. 0.0d0)  wall = wall + dble(cmax)/dble(crate)
c
c     get elapsed CPU time via Fortran intrinsic functions
c
      tcpu0 = tcpu
      call cpu_time (time)
      tcpu = dble(time)
      cpu = tcpu - tcpu0
      return
      end
