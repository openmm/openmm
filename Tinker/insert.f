c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine insert  --  insert atom into coordinates list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "insert" adds the specified atom to the Cartesian
c     coordinates list and shifts the remaining atoms
c
c
      subroutine insert (iatom)
      use sizes
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      implicit none
      integer i,j,iatom
c
c
c     increase by one the total number of atoms
c
      n = n + 1
c
c     shift the atom coordinates, types and connectivities
c
      do i = n, iatom+1, -1
         name(i) = name(i-1)
         x(i) = x(i-1)
         y(i) = y(i-1)
         z(i) = z(i-1)
         type(i) = type(i-1)
         n12(i) = n12(i-1)
         do j = 1, n12(i)
            i12(j,i) = i12(j,i-1)
         end do
      end do
c
c     put new atom at the origin with a big atom type number
c
      name(iatom) = 'NEW'
      x(iatom) = 0.0d0
      y(iatom) = 0.0d0
      z(iatom) = 0.0d0
      type(iatom) = maxtyp + 1
      n12(iatom) = 0
c
c     shift the connected atom lists to allow the insertion
c
      do i = 1, n
         do j = 1, n12(i)
            if (i12(j,i) .ge. iatom) then
               i12(j,i) = i12(j,i) + 1
            end if
         end do
      end do
c
c     write a message to describe the atom insertion
c
      if (debug) then
         write (iout,10)  iatom
   10    format (' INSERT  --  Inserting Atom Number :',i8)
      end if
      return
      end
