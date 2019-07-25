c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program molxyz  --  MDL MOL file to Cartesian coordinates  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "molxyz" takes as input a MDL MOL coordinates file,
c     converts to and then writes out Cartesian coordinates
c
c
      program molxyz
      use files
      use iounit
      use titles
      implicit none
      integer ixyz,freeunit
      character*240 xyzfile
c
c
c     get and read the MDL MOL format file
c
      call initial
      call getmol
      write (iout,10)  title(1:ltitle)
   10 format (/,' Title :  ',a)
c
c     write out the Cartesian coordinates file
c
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
      call prtxyz (ixyz)
      close (unit=ixyz)
c
c     perform any final tasks before program exit
c
      call final
      end
