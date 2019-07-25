c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module bndstr  --  bond stretches in the current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     nbond   total number of bond stretches in the system
c     ibnd    numbers of the atoms in each bond stretch
c     bk      bond stretch force constants (kcal/mole/Ang**2)
c     bl      ideal bond length values in Angstroms
c
c
      module bndstr
      implicit none
      integer nbond
      integer, allocatable :: ibnd(:,:)
      real*8, allocatable :: bk(:)
      real*8, allocatable :: bl(:)
      save
      end
