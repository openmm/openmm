c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program prtmol2  --  output of Sybyl MOL2 coordinate file  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtmol2" writes out a set of coordinates in Sybyl MOL2
c     format to an external disk file
c
c
      subroutine prtmol2 (isyb)
      use sizes
      use atomid
      use atoms
      use bndstr
      use couple
      use files
      use iounit
      use titles
      implicit none
      integer i,j,k,isyb
      integer thousand,hundred
      integer tens,ones
      logical opened
      character*1 digit(0:9)
      character*4 atmtyp,number
      character*7 atmnam
      character*240 sybylfile
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     open output unit if not already done
c
      inquire (unit=isyb,opened=opened)
      if (.not. opened) then
         sybylfile = filename(1:leng)//'.mol2'
         call version (sybylfile,'new')
         open (unit=isyb,file=sybylfile,status='new')
      end if
c
c     write the molecule record type indicator
c
      write (isyb,10)
   10 format ('@<TRIPOS>MOLECULE')
      if (ltitle .eq. 0) then
         write (isyb,20)
   20    format ('****')
      else
         write (isyb,30)  title(1:ltitle)
   30    format (a)
      end if
      write (isyb,40)  n,nbond,1
   40 format (3i8)
      write (isyb,50)
   50 format ('SMALL')
      write (isyb,60)
   60 format ('NO_CHARGES')
c
c     write the atom record type indicator
c
      write (isyb,70)
   70 format (/,'@<TRIPOS>ATOM')
      do i = 1, n
c
c     set the Sybyl atom_name for the atom
c
         thousand = i / 1000
         hundred = (i - 1000*thousand) / 100
         tens = (i - 1000*thousand - 100*hundred) / 10
         ones = i - 1000*thousand - 100*hundred - 10*tens
         number(1:1) = digit(thousand)
         number(2:2) = digit(hundred)
         number(3:3) = digit(tens)
         number(4:4) = digit(ones)
         if (number(1:1) .eq. '0')  number(1:1) = ' '
         if (number(2:2).eq.'0' .and. number(1:1).eq.' ') then
            number(2:2) = ' '
         end if
         if (number(3:3).eq.'0' .and. number(2:2).eq.' ') then
            number(3:3) = ' '
         end if
         atmnam = name(i)//number
         do j = 1, 6
            do while (atmnam(j:j) .eq. ' ')
               do k = j, 6
                  atmnam(k:k) = atmnam(k+1:k+1)
               end do
               atmnam(7:7) = '*'
            end do
         end do
         do j = 1, 7
            if (atmnam(j:j) .eq. '*')  atmnam(j:j) = ' '
         end do
c
c     set the Sybyl atom_type for the atom
c
         atmtyp = name(i)//' '
         if (atmtyp .eq. 'C  ') then
            if (n12(i) .eq. 4)  atmtyp = 'C.3 '
            if (n12(i) .eq. 3)  atmtyp = 'C.2 '
            if (n12(i) .eq. 2)  atmtyp = 'C.1 '
         else if (atmtyp .eq. 'N  ') then
            if (n12(i) .ge. 3)  atmtyp = 'N.3 '
            if (n12(i) .eq. 2)  atmtyp = 'N.2 '
            if (n12(i) .eq. 1)  atmtyp = 'N.1 '
         else if (atmtyp .eq. 'N+ ') then
            atmtyp = 'N.4 '
         else if (atmtyp .eq. 'O  ') then
            if (n12(i) .ge. 2)  atmtyp = 'O.3 '
            if (n12(i) .le. 1)  atmtyp = 'O.2 '
         else if (atmtyp .eq. 'O- ') then
            atmtyp = 'O.2 '
         else if (atmtyp .eq. 'S  ') then
            if (n12(i) .ge. 2)  atmtyp = 'S.3 '
            if (n12(i) .le. 1)  atmtyp = 'S.2 '
         else if (atmtyp .eq. 'P  ') then
            atmtyp = 'P.3 '
         else if (atmtyp .eq. 'Lp ') then
            atmtyp = 'LP  '
         end if
         write (isyb,80)  i,atmnam,x(i),y(i),z(i),atmtyp
   80    format (i8,3x,a7,2x,3f12.6,3x,a4)
      end do
c
c     write the bond record type indicator
c
      write (isyb,90)
   90 format (/,'@<TRIPOS>BOND')
      do i = 1, nbond
         write (isyb,100)  i,(ibnd(j,i),j=1,2),1
  100    format (4i8)
      end do
c
c     write the substructure record type indicator
c
      write (isyb,110)
  110 format (/,'@<TRIPOS>SUBSTRUCTURE')
      write (isyb,120)  1,'****',1
  120 format (i8,12x,a4,i8)
      if (.not. opened)  close (unit=isyb)
      return
      end
