c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module ctran -- charge transfer term forcefield parameters   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     nct       total number ct active sites in the system
c     apre    prefactor a parameters in charge transfer 
c     bexp    exponential b parameters in charge transfer
c     ct2scale  ct parameter by which 1-2 ct energy scaled
c     ct3scale  ct parameter by which 1-3 ct energy scaled
c     ct4scale  ct parameter by which 1-4 ct energy scaled
c     ct5scale  ct parameter by which 1-5 ct energy scaled
c     aprerule  apre combining rule
c     bexprule  bexp combining rule
c
      module ctran
      use sizes
      implicit none
      integer nct
      integer, allocatable :: ict(:)
      integer, allocatable :: jct(:)
      real*8 apre(maxtyp)
      real*8 bexp(maxtyp)
      real*8 ct2scale,ct3scale
      real*8 ct4scale,ct5scale
      character*10 aprerule,bexprule
      save
      end
