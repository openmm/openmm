c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module urey  --  Urey-Bradley interactions in structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nurey   total number of Urey-Bradley terms in the system
c     iury    numbers of the atoms in each Urey-Bradley interaction
c     uk      Urey-Bradley force constants (kcal/mole/Ang**2)
c     ul      ideal 1-3 distance values in Angstroms
c
c
      module urey
      implicit none
      integer nurey
      integer, allocatable :: iury(:,:)
      real*8, allocatable :: uk(:)
      real*8, allocatable :: ul(:)
      save
      end
