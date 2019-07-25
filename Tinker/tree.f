c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module tree  --  potential smoothing search tree levels  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxpss   maximum number of potential smoothing levels
c
c     nlevel   number of levels of potential smoothing used
c     etree    energy reference value at the top of the tree
c     ilevel   smoothing deformation value at each tree level
c
c
      module tree
      implicit none
      integer maxpss
      parameter (maxpss=500)
      integer nlevel
      real*8 etree
      real*8 ilevel(0:maxpss)
      save
      end
