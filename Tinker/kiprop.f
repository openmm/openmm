c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kiprop  --  improper dihedral forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxndi   maximum number of improper dihedral parameter entries
c
c     dcon     force constant parameters for improper dihedrals
c     tdi      ideal dihedral angle values for improper dihedrals
c     kdi      string of atom classes for improper dihedral angles
c
c
      module kiprop
      implicit none
      integer maxndi
      parameter (maxndi=500)
      real*8 dcon(maxndi)
      real*8 tdi(maxndi)
      character*16 kdi(maxndi)
      save
      end
