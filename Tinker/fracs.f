c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module fracs  --  distances to molecular center of mass  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xfrac   fractional coordinate along a-axis of center of mass
c     yfrac   fractional coordinate along b-axis of center of mass
c     zfrac   fractional coordinate along c-axis of center of mass
c
c
      module fracs
      implicit none
      real*8, allocatable :: xfrac(:)
      real*8, allocatable :: yfrac(:)
      real*8, allocatable :: zfrac(:)
      save
      end
