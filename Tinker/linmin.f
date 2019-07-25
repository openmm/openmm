c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module linmin  --  line search minimization parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     stpmin   minimum step length in current line search direction
c     stpmax   maximum step length in current line search direction
c     cappa    stringency of line search (0=tight < cappa < 1=loose)
c     slpmax   projected gradient above which stepsize is reduced
c     angmax   maximum angle between search direction and -gradient
c     intmax   maximum number of interpolations during line search
c
c
      module linmin
      implicit none
      integer intmax
      real*8 stpmin
      real*8 stpmax
      real*8 cappa
      real*8 slpmax
      real*8 angmax
      save
      end
