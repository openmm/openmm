c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getxyz  --  get Cartesian coordinate structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getxyz" asks for a Cartesian coordinate file name,
c     then reads in the coordinates file
c
c
      subroutine getxyz
      use inform
      use iounit
      use output
      implicit none
      integer ixyz
      integer freeunit
      logical exist
      character*240 xyzfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Cartesian Coordinate File Name :  ',$)
         read (input,20)  xyzfile
   20    format (a240)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      close (unit=ixyz)
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' GETXYZ  --  Cartesian Coordinates File',
     &              ' does not Contain Any Atoms')
         call fatal
      end if
      return
      end
