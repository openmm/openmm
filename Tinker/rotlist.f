c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotlist  --  find atoms on one side of a bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotlist" generates the minimum list of all the atoms lying
c     to one side of a pair of directly bonded atoms; optionally
c     finds the minimal list by choosing the side with fewer atoms
c
c
      subroutine rotlist (base,partner)
      use sizes
      use atoms
      use couple
      use iounit
      use molcul
      use rotbnd
      use zclose
      implicit none
      integer i,k,ia,ib,swap
      integer base,partner
      integer mark,test,nattach
      integer, allocatable :: list(:)
      logical bonded
c
c
c     initialize the number of atoms to one side of the bond
c
      nrot = 0
c
c     remove any bonds needed for intramolecular ring closures
c
      do i = 1, nadd
         ia = iadd(1,i)
         ib = iadd(2,i)
         if (molcule(ia) .eq. molcule(ib)) then
            do k = 1, n12(ia)
               if (i12(k,ia) .eq. ib)  i12(k,ia) = 0
            end do
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. ia)  i12(k,ib) = 0
            end do
         end if
      end do
c
c     add any links needed to make intermolecular connections
c
      do i = 1, ndel
         ia = idel(1,i)
         ib = idel(2,i)
         if (molcule(ia) .ne. molcule(ib)) then
            if (n12(ia).eq.maxval .or. n12(ib).eq.maxval) then
               write (iout,10)
   10          format (/,' ROTLIST  --  Maximum Valence Exceeded;',
     &                    ' Increase MAXVAL')
               call fatal
            end if
            n12(ia) = n12(ia) + 1
            i12(n12(ia),ia) = ib
            n12(ib) = n12(ib) + 1
            i12(n12(ib),ib) = ia
         end if
      end do
c
c     check to see if the two atoms are still directly bonded
c
      bonded = .false.
      do i = 1, n12(base)
         if (i12(i,base) .eq. partner)  bonded = .true.
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(rot))  allocate (rot(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(0:n))
c
c     make a list of atoms to one side of this pair of atoms,
c     taking note of any rings in which the atom pair resides
c
      if (bonded) then
         list(0) = 1
         do i = 1, n
            rot(i) = 0
         end do
   20    continue
         nrot = 0
         do i = 1, n
            list(i) = 0
         end do
         list(base) = 1
         list(partner) = 1
         nattach = n12(base)
         do i = 1, nattach
            test = i12(i,base)
            if (list(test) .eq. 0) then
               nrot = nrot + 1
               if (use_short .and. nrot.ge.n/2)  goto 30
               rot(nrot) = test
               list(test) = 1
            end if
         end do
         do i = 1, n
            mark = rot(i)
            if (mark .eq. 0)  goto 40
            nattach = n12(mark)
            if (nattach .gt. 1) then
               do k = 1, nattach
                  test = i12(k,mark)
                  if (list(test) .eq. 0) then
                     nrot = nrot + 1
                     if (use_short .and. nrot.ge.n/2)  goto 30
                     rot(nrot) = test
                     list(test) = 1
                  end if
               end do
            end if
         end do
c
c     the list contains over half the total number of atoms,
c     so reverse the base and partner atoms, then start over
c
   30    continue
         swap = base
         base = partner
         partner = swap
         do i = 1, nrot
            rot(i) = 0
         end do
         goto 20
      end if
   40 continue
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     remove links added to make intermolecular connections
c
      do i = 1, ndel
         ia = idel(1,i)
         ib = idel(2,i)
         if (molcule(ia) .ne. molcule(ib)) then
            n12(ia) = n12(ia) - 1
            n12(ib) = n12(ib) - 1
         end if
      end do
c
c     add any bonds required for intramolecular ring closures
c
      do i = 1, nadd
         ia = iadd(1,i)
         ib = iadd(2,i)
         if (molcule(ia) .eq. molcule(ib)) then
            do k = 1, n12(ia)
               if (i12(k,ia) .eq. 0) then
                  i12(k,ia) = ib
                  goto 50
               end if
            end do
   50       continue
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. 0) then
                  i12(k,ib) = ia
                  goto 60
               end if
            end do
   60       continue
         end if
      end do
      return
      end
