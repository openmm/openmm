c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erxnfld1  --  reaction field energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erxnfld1" calculates the macroscopic reaction field energy
c     and derivatives with respect to Cartesian coordinates
c
c
      subroutine erxnfld1
      use sizes
      use atoms
      use deriv
      use energi
      implicit none
      integer i,j
c
c
c     zero out macroscopic reaction field energy and derivatives
c
      er = 0.0d0
      do i = 1, n
         do j = 1, 3
            der(j,i) = 0.0d0
         end do
      end do
      return
      end
