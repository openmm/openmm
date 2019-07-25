c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module chrono  --  clock time values for current program  ##
c     ##                                                            ##
c     ################################################################
c
c
c     twall   current processor wall clock time in seconds
c     tcpu    elapsed cpu time from start of program in seconds
c
c
      module chrono
      implicit none
      real*8 twall
      real*8 tcpu
      save
      end
