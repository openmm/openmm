c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module scales  --  optimization parameter scale factors  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     scale      multiplicative factor for each optimization parameter
c     set_scale  logical flag to show if scale factors have been set
c
c
      module scales
      use sizes
      implicit none
      real*8 scale(3*maxatm)
      logical set_scale
      save
      end
