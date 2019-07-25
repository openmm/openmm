c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bonds  --  locate and store covalent bonds  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bonds" finds the total number of covalent bonds and
c     stores the atom numbers of the atoms defining each bond
c
c
      subroutine bonds
      use sizes
      use atmlst
      use atoms
      use bndstr
      use couple
      use iounit
      implicit none
      integer i,j,k,m
      integer maxbnd
c
c
c     perform dynamic allocation of some global arrays
c
      maxbnd = 4 * n
      if (allocated(ibnd))  deallocate (ibnd)
      if (allocated(bndlist))  deallocate (bndlist)
      allocate (ibnd(2,maxbnd))
      allocate (bndlist(maxval,n))
c
c     loop over all atoms, storing the atoms in each bond
c
      nbond = 0
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            if (i .lt. k) then
               nbond = nbond + 1
               if (nbond .gt. maxbnd) then
                  write (iout,10)
   10             format (/,' BONDS  --  Too many Bonds; Increase',
     &                       ' MAXBND')
                  call fatal
               end if
               ibnd(1,nbond) = i
               ibnd(2,nbond) = k
               bndlist(j,i) = nbond
               do m = 1, n12(k)
                  if (i .eq. i12(m,k)) then
                     bndlist(m,k) = nbond
                     goto 20
                  end if
               end do
   20          continue
            end if
         end do
      end do
      return
      end
