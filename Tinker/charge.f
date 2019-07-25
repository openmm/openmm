c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module charge  --  partial charges in current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nion      total number of partial charges in system
c     iion      number of the atom site for each partial charge
c     jion      neighbor generation site for each partial charge
c     kion      cutoff switching site for each partial charge
c     pchg      magnitude of the partial charges (e-)
c     penalpha  charge penetration parameters 
c
c
      module charge
      implicit none
      integer nion
      integer, allocatable :: iion(:)
      integer, allocatable :: jion(:)
      integer, allocatable :: kion(:)
      real*8, allocatable :: pchg(:)
      save
      end
