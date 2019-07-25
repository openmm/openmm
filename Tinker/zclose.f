c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module zclose  --  Z-matrix ring openings and closures  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nadd   number of added bonds between Z-matrix atoms
c     ndel   number of bonds between Z-matrix bonds to delete
c     iadd   numbers of the atom pairs defining added bonds
c     idel   numbers of the atom pairs defining deleted bonds
c
c
      module zclose
      use sizes
      implicit none
      integer nadd,ndel
      integer iadd(2,maxatm)
      integer idel(2,maxatm)
      save
      end
