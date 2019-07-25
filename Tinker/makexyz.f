c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine makexyz  --  convert internal to Cartesian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "makexyz" generates a complete set of Cartesian coordinates
c     for a full structure from the internal coordinate values
c
c
      subroutine makexyz
      use sizes
      use atoms
      use zcoord
      implicit none
      integer i,chiral
      integer ia,ib,ic
      real*8 bond
      real*8 angle1
      real*8 angle2
c
c
c     loop over each atom in turn, finding its coordinates
c
      do i = 1, n
         ia = iz(1,i)
         ib = iz(2,i)
         ic = iz(3,i)
         chiral = iz(4,i)
         bond = zbond(i)
         angle1 = zang(i)
         angle2 = ztors(i)
         call xyzatm (i,ia,bond,ib,angle1,ic,angle2,chiral)
      end do
      return
      end
