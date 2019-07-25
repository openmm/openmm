c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine chkring  --  check atom set for small rings  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "chkring" tests an atom or a set of connected atoms for
c     their presence within a single 3- to 6-membered ring
c
c
      subroutine chkring (iring,ia,ib,ic,id)
      use couple
      implicit none
      integer i,j,k,m,p,q,r
      integer ia,ib,ic,id
      integer iring,nset
c
c
c     initialize the ring size and number of atoms to test
c
      iring = 0
      nset = 0
      if (ia .gt. 0)  nset = 1
      if (ib .gt. 0)  nset = 2
      if (ic .gt. 0)  nset = 3
      if (id .gt. 0)  nset = 4
c
c     cannot be in a ring if the terminal atoms are univalent
c
      if (nset .eq. 1) then
         if (n12(ia) .le. 1)  nset = 0
      else if (nset .eq. 2) then
         if (min(n12(ia),n12(ib)) .le. 1)  nset = 0
      else if (nset .eq. 3) then
         if (min(n12(ia),n12(ic)) .le. 1)  nset = 0
      else if (nset .eq. 4) then
         if (min(n12(ia),n12(id)) .le. 1)  nset = 0
      end if
c
c     check the input atoms for sequential connectivity
c
      if (nset .gt. 1) then
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .eq. i) then
               if (nset .eq. 2)  goto 10
               do k = 1, n12(ib)
                  m = i12(k,ib)
                  if (ic .eq. m) then
                     if (nset .eq. 3)  goto 10
                     do p = 1, n12(ic)
                        q = i12(p,ic)
                        if (id .eq. q)  goto 10
                     end do
                  end if
               end do
            end if
         end do
         nset = 0
   10    continue
      end if
c
c     check for an atom contained inside a small ring
c
      if (nset .eq. 1) then
         do j = 1, n12(ia)-1
            i = i12(j,ia)
            do k = j+1, n12(ia)
               m = i12(k,ia)
               do p = 1, n12(i)
                  if (m .eq. i12(p,i)) then
                     iring = 3
                     goto 20
                  end if
               end do
            end do
         end do
         do j = 1, n12(ia)-1
            i = i12(j,ia)
            do k = j+1, n12(ia)
               m = i12(k,ia)
               do p = 1, n12(i)
                  r = i12(p,i)
                  if (r .ne. ia) then
                     do q = 1, n12(m)
                        if (r .eq. i12(q,m)) then
                           iring = 4
                           goto 20
                        end if
                     end do
                  end if
               end do
            end do
         end do
         do j = 1, n13(ia)-1
            i = i13(j,ia)
            do k = j+1, n13(ia)
               m = i13(k,ia)
               do p = 1, n12(i)
                  if (m .eq. i12(p,i)) then
                     iring = 5
                     goto 20
                  end if
               end do
               do p = 1, n13(i)
                  if (m .eq. i13(p,i)) then
                     iring = 6
                     goto 20
                  end if
               end do
            end do
         end do
   20    continue
c
c     check for a bond contained inside a small ring
c
      else if (nset .eq. 2) then
         do j = 1, n12(ia)
            i = i12(j,ia)
            do k = 1, n12(ib)
               if (i .eq. i12(k,ib)) then
                  iring = 3
                  goto 30
               end if
            end do
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n12(ib)
                  m = i12(k,ib)
                  if (ia .ne. m) then
                     do p = 1, n12(i)
                        if (m .eq. i12(p,i)) then
                           iring = 4
                           goto 30
                        end if
                     end do
                  end if
               end do
            end if
         end do
         do j = 1, n13(ia)
            i = i13(j,ia)
            do k = 1, n13(ib)
               if (i .eq. i13(k,ib)) then
                  iring = 5
                  goto 30
               end if
            end do
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n13(ib)
                  m = i13(k,ib)
                  do p = 1, n13(i)
                     if (m .eq. i13(p,i)) then
                        iring = 6
                        do q = 1, n12(ia)
                           if (m .eq. i12(q,ia))  iring = 0
                        end do
                        if (iring .eq. 6)  goto 30
                     end if
                  end do
               end do
            end if
         end do
   30    continue
c
c     check for an angle contained inside a small ring
c
      else if (nset .eq. 3) then
         do j = 1, n12(ia)
            if (ic .eq. i12(j,ia)) then
               iring = 3
               goto 40
            end if
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n12(ic)
                  if (i .eq. i12(k,ic)) then
                     iring = 4
                     goto 40
                  end if
               end do
            end if
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n13(ic)
                  if (i .eq. i13(k,ic)) then
                     iring = 5
                     goto 40
                  end if
               end do
            end if
         end do
         do j = 1, n13(ia)
            i = i13(j,ia)
            if (ic .ne. i) then
               do k = 1, n13(ic)
                  if (i .eq. i13(k,ic)) then
                     iring = 6
                     goto 40
                  end if
               end do
            end if
         end do
   40    continue
c
c     check for a torsion contained inside a small ring
c
      else if (nset .eq. 4) then
         do j = 1, n12(ia)
            if (id .eq. i12(j,ia)) then
               iring = 4
               goto 50
            end if
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n12(id)
                  if (i .eq. i12(k,id)) then
                     iring = 5
                     goto 50
                  end if
               end do
            end if
         end do
         do j = 1, n12(ia)
            i = i12(j,ia)
            if (ib .ne. i) then
               do k = 1, n13(id)
                  if (i .eq. i13(k,id)) then
                     iring = 6
                     goto 50
                  end if
               end do
            end if
         end do
   50    continue
      end if
      return
      end
