c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module molcul  --  individual molecules in current system  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nmol      total number of separate molecules in the system
c     imol      first and last atom of each molecule in the list
c     kmol      contiguous list of the atoms in each molecule
c     molcule   number of the molecule to which each atom belongs
c     totmass   total weight of all the molecules in the system
c     molmass   molecular weight for each molecule in the system
c
c
      module molcul
      implicit none
      integer nmol
      integer, allocatable :: imol(:,:)
      integer, allocatable :: kmol(:)
      integer, allocatable :: molcule(:)
      real*8 totmass
      real*8, allocatable :: molmass(:)
      save
      end
