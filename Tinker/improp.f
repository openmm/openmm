c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module improp  --  improper dihedrals in current structure  ##
c     ##                                                              ##
c     ######################################################3###########
c
c
c     niprop   total number of improper dihedral angles in the system
c     iiprop   numbers of the atoms in each improper dihedral angle
c     kprop    force constant values for improper dihedral angles
c     vprop    ideal improper dihedral angle value in degrees
c
c
      module improp
      implicit none
      integer niprop
      integer, allocatable :: iiprop(:,:)
      real*8, allocatable :: kprop(:)
      real*8, allocatable :: vprop(:)
      save
      end
