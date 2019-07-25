c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module pitors  --  pi-system torsions in current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     npitors   total number of pi-system torsional interactions
c     ipit      numbers of the atoms in each pi-system torsion
c     kpit      2-fold pi-system torsional force constants
c
c
      module pitors
      implicit none
      integer npitors
      integer, allocatable :: ipit(:,:)
      real*8, allocatable :: kpit(:)
      save
      end
