c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module usolve  --  induced dipole preconditioner inverse  ##
c     ##                                                            ##
c     ################################################################
c
c
c     mindex   index into preconditioner inverse for CG solver
c     minv     preconditioner inverse for induced dipole CG solver
c
c
      module usolve
      implicit none
      integer, allocatable :: mindex(:)
      real*8, allocatable :: minv(:)
      save
      end
