c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra1" calculates any additional user defined potential
c     energy contribution and its first derivatives
c
c
      subroutine extra1
      use sizes
      use atoms
      use deriv
      use energi
      implicit none
      integer i
c
c
c     zero out the extra energy term and first derivatives
c
      ex = 0.0d0
      do i = 1, n
         dex(1,i) = 0.0d0
         dex(2,i) = 0.0d0
         dex(3,i) = 0.0d0
      end do
c
c     add any user-defined extra potentials and derivatives;
c     also increment intermolecular energy and virial as needed
c
c     e = ......
c     ex = ex + e
c     do i = 1, n
c        dex(1,i) = ......
c        dex(2,i) = ......
c        dex(3,i) = ......
c     end do
c
      return
      end
