c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine getmol2  --  get a Sybyl MOL2 format file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "getmol2" asks for a Sybyl MOL2 molecule file name,
c     then reads the coordinates from the file
c
c
      subroutine getmol2
      use files
      use iounit
      implicit none
      integer isyb
      integer freeunit
      logical exist
      character*240 sybylfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (sybylfile,exist)
      if (exist) then
         call basefile (sybylfile)
         call suffix (sybylfile,'mol2','old')
         inquire (file=sybylfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter a Sybyl MOL2 File Name :  ',$)
         read (input,20)  sybylfile
   20    format (a240)
         call basefile (sybylfile)
         call suffix (sybylfile,'mol2','old')
         inquire (file=sybylfile,exist=exist)
      end do
c
c     first open and then read the Sybyl MOL2 coordinates file
c
      isyb = freeunit ()
      open (unit=isyb,file=sybylfile,status='old')
      rewind (unit=isyb)
      call readmol2 (isyb)
      close (unit=isyb)
      return
      end
