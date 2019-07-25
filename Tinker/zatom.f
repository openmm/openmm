c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine zatom  --  adds a single atom to Z-matrix  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "zatom" adds an atom to the end of the current Z-matrix
c     and then increments the atom counter; atom type, defining
c     atoms and internal coordinates are passed as arguments
c
c
      subroutine zatom (bionum,bond,angle,dihed,iz1,iz2,iz3,iz4)
      use sizes
      use angbnd
      use atomid
      use atoms
      use bndstr
      use fields
      use iounit
      use kangs
      use katoms
      use kbonds
      use zclose
      use zcoord
      implicit none
      integer i,size
      integer bionum
      integer na,nb
      integer iz1,iz2
      integer iz3,iz4
      integer ita,itb
      integer itc,itd
      real*8 bond,angle
      real*8 dihed
      logical lookup
      character*4 pa,pb,pc,pd
      character*8 blank8,ptb
      character*12 blank12,pta
c
c
c     choose ideal or force field bond and angle values
c
      lookup = .false.
c
c     fill various arrays with information for this atom
c
      if (bionum .gt. 0) then
         type(n) = biotyp(bionum)
         if (type(n) .gt. 0) then
            name(n) = symbol(type(n))
         else
            name(n) = 'XXX'
         end if
         zbond(n) = bond
         zang(n) = angle
         ztors(n) = dihed
         if (ztors(n) .lt. -180.0d0) then
            ztors(n) = ztors(n) + 360.0d0
         else if (ztors(n) .gt. 180.0d0) then
            ztors(n) = ztors(n) - 360.0d0
         end if
         iz(1,n) = iz1
         iz(2,n) = iz2
         iz(3,n) = iz3
         iz(4,n) = iz4
c
c     find force field bond length and angle values
c
         if (lookup) then
            ita = 0
            itb = 0
            itc = 0
            itd = 0
            if (n .ne. 0)  ita = atmcls(type(n))
            if (iz1 .ne. 0)  itb = atmcls(type(iz1))
            if (iz2 .ne. 0)  itc = atmcls(type(iz2))
            if (iz3 .ne. 0)  itd = atmcls(type(iz3))
            blank8 = '        '
            blank12 = '            '
            do i = maxnb, 1, -1
               if (kb(i) .eq. blank8)  nb = i - 1
            end do
            do i = maxna, 1, -1
               if (ka(i) .eq. blank12)  na = i - 1
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               ptb = pa//pb
            else
               ptb = pb//pa
            end if
            do i = 1, nb
               if (kb(i) .eq. ptb) then
                  if (blen(i) .ne. 0.0d0)  zbond(n) = blen(i)
                  goto 10
               end if
            end do
   10       continue
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pta = pa//pb//pc
            else
               pta = pc//pb//pa
            end if
            do i = 1, na
               if (ka(i) .eq. pta) then
                  if (ang(1,i) .ne. 0.0d0)  zang(n) = ang(1,i)
                  goto 20
               end if
            end do
   20       continue
            if (iz4 .ne. 0) then
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itd,pd,size)
               if (ita .le. itd) then
                  pta = pa//pb//pd
               else
                  pta = pd//pb//pa
               end if
               do i = 1, na
                  if (ka(i) .eq. pta) then
                     if (ang(1,i) .ne. 0.0d0)  ztors(n) = ang(1,i)
                     goto 30
                  end if
               end do
   30          continue
            end if
         end if
c
c     increment atom counter and check for too many atoms
c
         n = n + 1
         if (n .gt. maxatm) then
            write (iout,40)  maxatm
   40       format (/,' ZATOM  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     add an extra bond to make a ring closure
c
      else if (bionum .eq. -1) then
         nadd = nadd + 1
         iadd(1,nadd) = iz1
         iadd(2,nadd) = iz2
c
c     delete an extra bond to make separate molecules
c
      else if (bionum .eq. -2) then
         ndel = ndel + 1
         idel(1,ndel) = iz1
         idel(2,ndel) = iz2
      end if
      return
      end
