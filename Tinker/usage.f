c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module usage  --  atoms active during energy computation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nuse   total number of active atoms in energy calculation
c     iuse   numbers of the atoms active in energy calculation
c     use    true if an atom is active, false if inactive
c
c
      module usage
      implicit none
      integer nuse
      integer, allocatable :: iuse(:)
      logical, allocatable :: use(:)
      save
      end
