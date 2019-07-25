c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program xyzsybyl  --  convert Cartesian to Sybyl MOL2 file  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "xyzsybyl" takes as input a Cartesian coordinates file,
c     converts to and then writes out a Sybyl MOL2 file
c
c
      program xyzsybyl
      use sizes
      use files
      use iounit
      use titles
      implicit none
      integer isyb,freeunit
      character*240 sybylfile
c
c
c     get and read the Cartesian coordinates file
c
      call initial
      call getxyz
      write (iout,10)  title(1:ltitle)
   10 format (' Title :  ',a)
c
c     get a list of the bonds in the structure
c
      call bonds
c
c     open a new version of the Sybyl MOL2 file
c
      isyb = freeunit ()
      sybylfile = filename(1:leng)//'.mol2'
      call version (sybylfile,'new')
      open (unit=isyb,file=sybylfile,status='new')
c
c     output the coordinates into Sybyl MOL2 format
c
      call prtmol2 (isyb)
      close (unit=isyb)
c
c     perform any final tasks before program exit
c
      call final
      end
