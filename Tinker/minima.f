c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module minima  --  general parameters for minimizations  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     fctmin    value below which function is deemed optimized
c     hguess    initial value for the H-matrix diagonal elements
c     maxiter   maximum number of iterations during optimization
c     nextiter  iteration number to use for the first iteration
c
c
      module minima
      implicit none
      integer maxiter
      integer nextiter
      real*8 fctmin
      real*8 hguess
      save
      end
