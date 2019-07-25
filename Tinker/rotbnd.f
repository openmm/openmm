c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module rotbnd  --  molecule partitions for bond rotation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nrot        total number of atoms moving when bond rotates
c     rot         atom numbers of atoms moving when bond rotates
c     use_short   logical flag governing use of shortest atom list
c
c
      module rotbnd
      implicit none
      integer nrot
      integer, allocatable :: rot(:)
      logical use_short
      save
      end
