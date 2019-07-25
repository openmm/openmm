c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine connect  --  attached atom list from Z-matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "connect" sets up the attached atom arrays
c     starting from a set of internal coordinates
c
c
      subroutine connect
      use sizes
      use atoms
      use couple
      use zcoord
      use zclose
      implicit none
      integer i,j,k
      integer id1,id2
c
c
c     zero out the number of atoms attached to each atom
c
      do i = 1, n
         n12(i) = 0
      end do
c
c     loop over the bonds in the Z-matrix, adding each bond
c     to the attach atom lists unless it is to be removed
c
      do i = 2, n
         k = iz(1,i)
         do j = 1, ndel
            id1 = idel(1,j)
            id2 = idel(2,j)
            if ((i.eq.id1 .and. k.eq.id2) .or.
     &          (i.eq.id2 .and. k.eq.id1))  goto 10
         end do
         n12(i) = n12(i) + 1
         n12(k) = n12(k) + 1
         i12(n12(i),i) = k
         i12(n12(k),k) = i
   10    continue
      end do
c
c     add any extra bonds used to make ring closures
c
      do i = 1, nadd
         do j = 1, 2
            k = iadd(j,i)
            n12(k) = n12(k) + 1
            i12(n12(k),k) = iadd(3-j,i)
         end do
      end do
c
c     sort the attached atom lists into ascending order
c
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
      return
      end
