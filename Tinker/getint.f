c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getint  --  get internal coordinate structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getint" asks for an internal coordinate file name, then reads
c     the internal coordinates and computes Cartesian coordinates
c
c
      subroutine getint
      use sizes
      use atoms
      use inform
      use iounit
      use output
      implicit none
      integer izmt
      integer freeunit
      logical exist
      logical clash
      character*240 intfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (intfile,exist)
      if (exist) then
         call basefile (intfile)
         call suffix (intfile,'int','old')
         inquire (file=intfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Internal Coordinate File Name :  ',$)
         read (input,20)  intfile
   20    format (a240)
         call basefile (intfile)
         call suffix (intfile,'int','old')
         inquire (file=intfile,exist=exist)
      end do
c
c     first open and then read the internal coordinates file
c
      coordtype = 'INTERNAL'
      izmt = freeunit ()
      open (unit=izmt,file=intfile,status='old')
      rewind (unit=izmt)
      call readint (izmt)
      close (unit=izmt)
c
c     quit if the internal coordinates file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' GETINT  --  Internal Coordinates File',
     &              ' does not Contain Any Atoms')
         call fatal
      end if
c
c     convert internal to Cartesian coordinates
c
      call connect
      call makexyz
c
c     check for atoms with identical coordinates
c
      clash = .false.
      if (n .le. 10000)  call chkxyz (clash)
      return
      end
