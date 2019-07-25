c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kanang  --  angle-angle term forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     anan   angle-angle cross term parameters for each atom class
c
c
      module kanang
      use sizes
      implicit none
      real*8 anan(3,maxclass)
      save
      end
