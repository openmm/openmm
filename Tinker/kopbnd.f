c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kopbnd  --  out-of-plane bend forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnopb   maximum number of out-of-plane bending entries
c
c     opbn      force constant parameters for out-of-plane bending
c     kopb      string of atom classes for out-of-plane bending
c
c
      module kopbnd
      implicit none
      integer maxnopb
      parameter (maxnopb=500)
      real*8 opbn(maxnopb)
      character*16 kopb(maxnopb)
      save
      end
