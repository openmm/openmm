c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module angtor  --  angle-torsions in current structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nangtor   total number of angle-torsion interactions
c     iat       torsion and angle numbers used in angle-torsion
c     kant      1-, 2- and 3-fold angle-torsion force constants
c
c
      module angtor
      implicit none
      integer nangtor
      integer, allocatable :: iat(:,:)
      real*8, allocatable :: kant(:,:)
      save
      end
