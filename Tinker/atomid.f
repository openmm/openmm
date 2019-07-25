c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module atomid  --  atomic properties for current atoms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tag       integer atom labels from input coordinates file
c     class     atom class number for each atom in the system
c     atomic    atomic number for each atom in the system
c     valence   valence number for each atom in the system
c     mass      atomic weight for each atom in the system
c     name      atom name for each atom in the system
c     story     descriptive type for each atom in system
c
c
      module atomid
      use sizes
      implicit none
      integer tag(maxatm)
      integer class(maxatm)
      integer atomic(maxatm)
      integer valence(maxatm)
      real*8 mass(maxatm)
      character*3 name(maxatm)
      character*24 story(maxatm)
      save
      end
