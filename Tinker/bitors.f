c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine bitors  --  locate and store bitorsions  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "bitors" finds the total number of bitorsions as pairs
c     of adjacent torsional angles, and the numbers of the five
c     atoms defining each bitorsion
c
c
      subroutine bitors
      use sizes
      use angbnd
      use atoms
      use bitor
      use couple
      use iounit
      implicit none
      integer i,j,k
      integer ia,ib,ic,id,ie
      integer maxbitor
c
c
c     perform dynamic allocation of some global arrays
c
      maxbitor = 54 * n
      if (allocated(ibitor))  deallocate (ibitor)
      allocate (ibitor(5,maxbitor))
c
c     loop over all angles, storing the atoms in each bitorsion
c
      nbitor = 0
      do i = 1, nangle
         ib = iang(1,i)
         ic = iang(2,i)
         id = iang(3,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id) then
               do k = 1, n12(id)
                  ie = i12(k,id)
                  if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                     nbitor = nbitor + 1
                     if (nbitor .gt. maxbitor) then
                        write (iout,10)
   10                   format (/,' BITORS  --  Too many Adjacent',
     &                             ' Torsions; Increase MAXBITOR')
                        call fatal
                     end if
                     ibitor(1,nbitor) = ia
                     ibitor(2,nbitor) = ib
                     ibitor(3,nbitor) = ic
                     ibitor(4,nbitor) = id
                     ibitor(5,nbitor) = ie
                  end if
               end do
            end if
         end do
      end do
      return
      end
