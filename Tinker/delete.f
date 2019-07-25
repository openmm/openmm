c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine delete  --  remove atom from coordinates list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "delete" removes a specified atom from the Cartesian
c     coordinates list and shifts the remaining atoms
c
c
      subroutine delete (iatom)
      use sizes
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      implicit none
      integer i,j,k,m,iatom
c
c
c     reduce by one the total number of atoms
c
      n = n - 1
c
c     shift the atom coordinates, types and connectivities
c
      do i = iatom, n
         name(i) = name(i+1)
         x(i) = x(i+1)
         y(i) = y(i+1)
         z(i) = z(i+1)
         type(i) = type(i+1)
         n12(i) = n12(i+1)
         do j = 1, n12(i)
            i12(j,i) = i12(j,i+1)
         end do
      end do
c
c     remove connections to deleted atom and shift the lists
c
      do i = 1, n
         m = 0
         do j = 1, n12(i)
            if (i12(j,i) .eq. iatom) then
               m = m + 1
               do k = j, n12(i)-1
                  i12(k,i) = i12(k+1,i)
               end do
            end if
         end do
         n12(i) = n12(i) - m
         do j = 1, n12(i)
            if (i12(j,i) .gt. iatom)  i12(j,i) = i12(j,i) - 1
         end do
      end do
c
c     write a message to describe the atom deletion
c
      if (debug) then
         write (iout,10)  iatom
   10    format (' DELETE  --  Deleting Atom Number :',i8)
      end if
      return
      end
