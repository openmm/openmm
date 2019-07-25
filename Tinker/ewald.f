c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ewald  --  Ewald summation parameters and options  ##
c     ##                                                            ##
c     ################################################################
c
c
c     aewald     Ewald convergence coefficient value (Ang-1)
c     boundary   Ewald boundary condition; none, tinfoil or vacuum
c
c
      module ewald
      implicit none
      real*8 aewald
      character*7 boundary
      save
      end
