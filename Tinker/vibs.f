c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module vibs  --  iterative vibrational analysis components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     phi        trial vectors for iterative vibrational analysis
c     phik       alternate vectors for iterative vibrational analysis
c     pwork      temporary work array for eigenvector transformation
c
c
      module vibs
      implicit none
      real*8, allocatable :: phi(:,:)
      real*8, allocatable :: phik(:,:)
      real*8, allocatable :: pwork(:,:)
      save
      end
