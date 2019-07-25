c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module ring  --  number and location of ring structures  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nring3   total number of 3-membered rings in the system
c     nring4   total number of 4-membered rings in the system
c     nring5   total number of 5-membered rings in the system
c     nring6   total number of 6-membered rings in the system
c     iring3   numbers of the atoms involved in each 3-ring
c     iring4   numbers of the atoms involved in each 4-ring
c     iring5   numbers of the atoms involved in each 5-ring
c     iring6   numbers of the atoms involved in each 6-ring
c
c
      module ring
      implicit none
      integer nring3
      integer nring4
      integer nring5
      integer nring6
      integer, allocatable :: iring3(:,:)
      integer, allocatable :: iring4(:,:)
      integer, allocatable :: iring5(:,:)
      integer, allocatable :: iring6(:,:)
      save
      end
