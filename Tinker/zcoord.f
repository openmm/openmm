c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module zcoord  --  Z-matrix internal coordinate values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     iz      defining atom numbers for each Z-matrix atom
c     zbond   bond length used to define each Z-matrix atom
c     zang    bond angle used to define each Z-matrix atom
c     ztors   angle or torsion used to define Z-matrix atom
c
c
      module zcoord
      use sizes
      implicit none
      integer iz(4,maxatm)
      real*8 zbond(maxatm)
      real*8 zang(maxatm)
      real*8 ztors(maxatm)
      save
      end
