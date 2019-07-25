c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module tortor  --  torsion-torsions in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     ntortor   total number of torsion-torsion interactions
c     itt       atoms and parameter indices for torsion-torsion
c
c
      module tortor
      implicit none
      integer ntortor
      integer, allocatable :: itt(:,:)
      save
      end
