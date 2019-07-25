c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readint  --  input of internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readint" gets a set of Z-matrix internal coordinates
c     from an external file
c
c
      subroutine readint (izmt)
      use sizes
      use atomid
      use atoms
      use files
      use inform
      use iounit
      use titles
      use zclose
      use zcoord
      implicit none
      integer i,j,izmt
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      logical exist,opened
      logical quit
      character*240 intfile
      character*240 record
      character*240 string
c
c
c     initialize the total number of atoms in the system
c
      n = 0
c
c     open the input file if it has not already been done
c
      inquire (unit=izmt,opened=opened)
      if (.not. opened) then
         intfile = filename(1:leng)//'.int'
         open (unit=izmt,file=intfile,status='old')
         rewind (unit=izmt)
         call version (intfile,'old')
         inquire (file=intfile,exist=exist)
         if (exist) then
            open (unit=izmt,file=intfile,status='old')
            rewind (unit=izmt)
         else
            write (iout,10)
   10       format (/,' READINT  --  Unable to Find the Internal',
     &                 ' Coordinates File')
            call fatal
         end if
      end if
c
c     read first line and return if already at end of file
c
      quit = .false.
      abort = .true.
      size = 0
      do while (size .eq. 0)
         read (izmt,20,err=70,end=70)  record
   20    format (a240)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=70,end=70)  n
c
c     extract the title and determine its length
c
      string = record(next:240)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         title = ' '
         ltitle = 0
      else
         title = string(first:last)
         ltitle = trimtext (title)
      end if
c
c     check for too many total atoms in the file
c
      if (n .le. 0) then
         write (iout,30)
   30    format (/,' READINT  --  The Coordinate File Does Not',
     &              ' Contain Any Atoms')
         call fatal
      else if (n .gt. maxatm) then
         write (iout,40)  maxatm
   40    format (' READINT  --  The Maximum of',i9,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = 0
         name(i) = '   '
         zbond(i) = 0.0d0
         zang(i) = 0.0d0
         ztors(i) = 0.0d0
         type(i) = 0
         do j = 1, 4
            iz(j,i) = 0
         end do
      end do
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         size = 0
         do while (size .eq. 0)
            read (izmt,50,err=70,end=70)  record
   50       format (a240)
            size = trimtext (record)
         end do
         read (record,*,err=70,end=70)  tag(i)
         call getword (record,name(i),next)
         string = record(next:240)
         read (string,*,err=60,end=60)  type(i),iz(1,i),zbond(i),
     &                                  iz(2,i),zang(i),iz(3,i),
     &                                  ztors(i),iz(4,i)
   60    continue
      end do
      quit = .false.
   70 continue
      if (.not. opened)  close (unit=izmt)
c
c     an error occurred in reading the Z-matrix coordinates
c
      if (quit) then
         write (iout,80)  i
   80    format (' READZ  --  Error in Z-Matrix File at Atom',i9)
         call fatal
      end if
c
c     read in any additional bonds to be added or deleted
c
      nadd = 0
      ndel = 0
      read (izmt,90,err=130,end=130)
   90 format ()
      do i = 1, maxatm
         read (izmt,100,err=130,end=130)  record
  100    format (a240)
         read (record,*,err=110,end=110)  (iadd(j,i),j=1,2)
         nadd = i
      end do
  110 continue
      do i = 1, maxatm
         read (izmt,120,err=130,end=130)  record
  120    format (a240)
         read (record,*,err=130,end=130)  (idel(j,i),j=1,2)
         ndel = i
      end do
  130 continue
      return
      end
