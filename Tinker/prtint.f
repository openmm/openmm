c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prtint  --  output of internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtint" writes out a set of Z-matrix internal
c     coordinates to an external disk file
c
c
      subroutine prtint (izmt)
      use sizes
      use atomid
      use atoms
      use files
      use inform
      use titles
      use zclose
      use zcoord
      implicit none
      integer i,k,izmt
      logical opened
      character*2 atmc
      character*5 bndc,angc
      character*43 fstr
      character*240 zmtfile
c
c
c     open the output unit if not already done
c
      inquire (unit=izmt,opened=opened)
      if (.not. opened) then
         zmtfile = filename(1:leng)//'.int'
         call version (zmtfile,'new')
         open (unit=izmt,file=zmtfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      if (digits .le. 6) then
         bndc = 'f10.5'
         angc = 'f10.4'
      else if (digits .le. 8) then
         bndc = 'f12.7'
         angc = 'f12.6'
      else
         bndc = 'f14.9'
         angc = 'f14.8'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (izmt,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (izmt,fstr(1:9))  n,title(1:ltitle)
      end if
c
c     output of first three atoms is handled separately
c
      fstr = '('//atmc//',2x,a3,i6,'//atmc//','//bndc//','//atmc//
     &          ','//angc//','//atmc//','//angc//','//'i6)'
      if (n .ge. 1)
     &   write (izmt,fstr)  1,name(1),type(1)
      if (n .ge. 2)
     &   write (izmt,fstr)  2,name(2),type(2),iz(1,2),zbond(2)
      if (n .ge. 3)
     &   write (izmt,fstr)  3,name(3),type(3),iz(1,3),zbond(3),
     &                      iz(2,3),zang(3)
c
c     convert torsional angles to lie in standard range
c
      do i = 4, n
         if (iz(4,i) .eq. 0) then
            do while (ztors(i) .lt. -180.0d0)
               ztors(i) = ztors(i) + 360.0d0
            end do
            do while (ztors(i) .gt. 180.0d0)
               ztors(i) = ztors(i) - 360.0d0
            end do
         end if
      end do
c
c     output the fourth through final atoms
c
      do i = 4, n
         write (izmt,fstr)  i,name(i),type(i),iz(1,i),zbond(i),
     &                      iz(2,i),zang(i),iz(3,i),ztors(i),iz(4,i)
      end do
c
c     addition and deletion of bonds as required
c
      if (nadd.ne.0 .or. ndel.ne.0) then
         fstr = '(2'//atmc//')'
         write (izmt,'()')
         do i = 1, nadd
            write (izmt,fstr(1:5))  (iadd(k,i),k=1,2)
         end do
         if (ndel .ne. 0)  write (izmt,'()')
         do i = 1, ndel
            write (izmt,fstr(1:5))  (idel(k,i),k=1,2)
         end do
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=izmt)
      return
      end
