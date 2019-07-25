c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module angang  --  angle-angles in current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     nangang   total number of angle-angle interactions
c     iaa       angle numbers used in each angle-angle term
c     kaa       force constant for angle-angle cross terms
c
c
      module angang
      implicit none
      integer nangang
      integer, allocatable :: iaa(:,:)
      real*8, allocatable :: kaa(:)
      save
      end
