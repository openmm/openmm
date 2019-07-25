c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "extra" calculates any additional user defined potential
c     energy contribution
c
c
      subroutine extra
      use energi
      implicit none
c
c
c     zero out the energy due to extra potential terms
c
      ex = 0.0d0
c
c     add any user-defined extra potentials below here
c
c     e = ......
c     ex = ex + e
c
      return
      end
