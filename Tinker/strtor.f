c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module strtor  --  stretch-torsions in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nstrtor   total number of stretch-torsion interactions
c     ist       torsion and bond numbers used in stretch-torsion
c     kst       1-, 2- and 3-fold stretch-torsion force constants
c
c
      module strtor
      implicit none
      integer nstrtor
      integer, allocatable :: ist(:,:)
      real*8, allocatable :: kst(:,:)
      save
      end
