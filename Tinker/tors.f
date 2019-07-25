c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module tors  --  torsional angles in current structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     ntors   total number of torsional angles in the system
c     itors   numbers of the atoms in each torsional angle
c     tors1   1-fold amplitude and phase for each torsional angle
c     tors2   2-fold amplitude and phase for each torsional angle
c     tors3   3-fold amplitude and phase for each torsional angle
c     tors4   4-fold amplitude and phase for each torsional angle
c     tors5   5-fold amplitude and phase for each torsional angle
c     tors6   6-fold amplitude and phase for each torsional angle
c
c
      module tors
      implicit none
      integer ntors
      integer, allocatable :: itors(:,:)
      real*8, allocatable :: tors1(:,:)
      real*8, allocatable :: tors2(:,:)
      real*8, allocatable :: tors3(:,:)
      real*8, allocatable :: tors4(:,:)
      real*8, allocatable :: tors5(:,:)
      real*8, allocatable :: tors6(:,:)
      save
      end
