c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kchrge  --  partial charge forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     chg   partial charge parameters for each atom type
c
c
      module kchrge
      implicit none
      real*8, allocatable :: chg(:)
      save
      end
