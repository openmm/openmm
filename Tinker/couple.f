c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module couple  --  atom neighbor connectivity lists   ##
c     ##                                                        ##
c     ############################################################
c
c
c     n12      number of atoms directly bonded to each atom
c     n13      number of atoms in a 1-3 relation to each atom
c     n14      number of atoms in a 1-4 relation to each atom
c     n15      number of atoms in a 1-5 relation to each atom
c     i12      atom numbers of atoms 1-2 connected to each atom
c     i13      atom numbers of atoms 1-3 connected to each atom
c     i14      atom numbers of atoms 1-4 connected to each atom
c     i15      atom numbers of atoms 1-5 connected to each atom
c
c
      module couple
      use sizes
      implicit none
      integer n12(maxatm)
      integer, allocatable :: n13(:)
      integer, allocatable :: n14(:)
      integer, allocatable :: n15(:)
      integer i12(maxval,maxatm)
      integer, allocatable :: i13(:,:)
      integer, allocatable :: i14(:,:)
      integer, allocatable :: i15(:,:)
      save
      end
