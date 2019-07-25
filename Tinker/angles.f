c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine angles  --  locate and store bond angles  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "angles" finds the total number of bond angles and stores
c     the atom numbers of the atoms defining each angle; for
c     each angle to a trivalent central atom, the third bonded
c     atom is stored for use in out-of-plane bending
c
c
      subroutine angles
      use sizes
      use angbnd
      use atmlst
      use atoms
      use couple
      use iounit
      implicit none
      integer i,j,k,m
      integer maxang
c
c
c     perform dynamic allocation of some global arrays
c
      maxang = 6 * n
      if (allocated(iang))  deallocate (iang)
      if (allocated(anglist))  deallocate (anglist)
      allocate (iang(4,maxang))
      allocate (anglist(maxval*(maxval-1)/2,n))
c
c     loop over all atoms, storing the atoms in each bond angle
c
      nangle = 0
      do i = 1, n
         m = 0
         do j = 1, n12(i)-1
            do k = j+1, n12(i)
               nangle = nangle + 1
               if (nangle .gt. maxang) then
                  write (iout,10)
   10             format (/,' ANGLES  --  Too many Bond Angles;',
     &                       ' Increase MAXANG')
                  call fatal
               end if
               m = m + 1
               anglist(m,i) = nangle
               iang(1,nangle) = i12(j,i)
               iang(2,nangle) = i
               iang(3,nangle) = i12(k,i)
               iang(4,nangle) = 0
            end do
         end do
c
c     set the out-of-plane atom for angles at trivalent centers
c
         if (n12(i) .eq. 3) then
            iang(4,nangle) = i12(1,i)
            iang(4,nangle-1) = i12(2,i)
            iang(4,nangle-2) = i12(3,i)
         end if
      end do
      return
      end
