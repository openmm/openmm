c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module atmlst  --  local geometry indices for each atom  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     bndlist   list of the bond numbers involving each atom
c     anglist   list of the angle numbers centered on each atom
c
c
      module atmlst
      implicit none
      integer, allocatable :: bndlist(:,:)
      integer, allocatable :: anglist(:,:)
      save
      end
