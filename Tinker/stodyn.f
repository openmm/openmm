c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module stodyn  --  SD trajectory frictional coefficients  ##
c     ##                                                            ##
c     ################################################################
c
c
c     friction    global frictional coefficient for exposed particle
c     fgamma      atomic frictional coefficients for each atom
c     use_sdarea  logical flag to use surface area friction scaling
c
c
      module stodyn
      implicit none
      real*8 friction
      real*8, allocatable :: fgamma(:)
      logical use_sdarea
      save
      end
