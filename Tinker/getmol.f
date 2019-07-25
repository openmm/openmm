c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine getmol  --  get a MDL MOL format file  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "getmol" asks for a MDL MOL molecule file name,
c     then reads the coordinates from the file
c
c
      subroutine getmol
      use files
      use iounit
      implicit none
      integer imdl
      integer freeunit
      logical exist
      character*240 mdlfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (mdlfile,exist)
      if (exist) then
         call basefile (mdlfile)
         call suffix (mdlfile,'mol','old')
         inquire (file=mdlfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter a MDL MOL File Name :  ',$)
         read (input,20)  mdlfile
   20    format (a240)
         call basefile (mdlfile)
         call suffix (mdlfile,'mol','old')
         inquire (file=mdlfile,exist=exist)
      end do
c
c     first open and then read the MDL MOL coordinates file
c
      imdl = freeunit ()
      open (unit=imdl,file=mdlfile,status='old')
      rewind (unit=imdl)
      call readmol (imdl)
      close (unit=imdl)
      return
      end
