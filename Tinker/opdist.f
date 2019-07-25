c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module opdist  --  out-of-plane distances in structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nopdist   total number of out-of-plane distances in the system
c     iopb      numbers of the atoms in each out-of-plane distance
c     opdk      force constant values for out-of-plane distance
c
c
      module opdist
      implicit none
      integer nopdist
      integer, allocatable :: iopd(:,:)
      real*8, allocatable :: opdk(:)
      save
      end
