c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module bitor  --  bitorsions in the current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nbitor  total number of bitorsions in the system
c     ibitor  numbers of the atoms in each bitorsion
c
c
      module bitor
      implicit none
      integer nbitor
      integer, allocatable :: ibitor(:,:)
      save
      end
