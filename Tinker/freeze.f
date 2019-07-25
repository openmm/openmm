c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module freeze  --  definition of holonomic constraints  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nrat         number of holonomic distance constraints to apply
c     nratx        number of atom group holonomic constraints to apply
c     iratx        group number of group in a holonomic constraint
c     kratx        spatial constraint type (1=plane, 2=line, 3=point)
c     irat         atom numbers of atoms in a holonomic constraint
c     rateps       convergence tolerance for holonomic constraints
c     krat         ideal distance value for holonomic constraint
c     use_rattle   logical flag to set use of holonomic contraints
c     ratimage     flag to use minimum image for holonomic constraint
c
c
      module freeze
      implicit none
      integer nrat,nratx
      integer, allocatable :: iratx(:)
      integer, allocatable :: kratx(:)
      integer, allocatable :: irat(:,:)
      real*8 rateps
      real*8, allocatable :: krat(:)
      logical use_rattle
      logical, allocatable :: ratimage(:)
      save
      end
