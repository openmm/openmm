c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module group  --  partitioning of system into atom groups  ##
c     ##                                                             ##
c     #################################################################
c
c
c     ngrp        total number of atom groups in the system
c     kgrp        contiguous list of the atoms in each group
c     grplist     number of the group to which each atom belongs
c     igrp        first and last atom of each group in the list
c     grpmass     total mass of all the atoms in each group
c     wgrp        weight for each set of group-group interactions
c     use_group   flag to use partitioning of system into groups
c     use_intra   flag to include only intragroup interactions
c     use_inter   flag to include only intergroup interactions
c
c
      module group
      implicit none
      integer ngrp
      integer, allocatable :: kgrp(:)
      integer, allocatable :: grplist(:)
      integer, allocatable :: igrp(:,:)
      real*8, allocatable :: grpmass(:)
      real*8, allocatable :: wgrp(:,:)
      logical use_group
      logical use_intra
      logical use_inter
      save
      end
