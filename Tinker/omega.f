c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module omega  --  torsional space dihedral angle values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nomega   number of dihedral angles allowed to rotate
c     iomega   numbers of two atoms defining rotation axis
c     zline    line number in Z-matrix of each dihedral angle
c     dihed    current value in radians of each dihedral angle
c
c
      module omega
      implicit none
      integer nomega
      integer, allocatable :: iomega(:,:)
      integer, allocatable :: zline(:)
      real*8, allocatable :: dihed(:)
      save
      end
