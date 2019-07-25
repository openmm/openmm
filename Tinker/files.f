c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module files  --  name & number of current structure file  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nprior     number of previously existing cycle files
c     ldir       length in characters of the directory name
c     leng       length in characters of the base filename
c     filename   base filename used by default for all files
c     outfile    output filename used for intermediate results
c
c
      module files
      implicit none
      integer nprior
      integer ldir,leng
      character*240 filename
      character*240 outfile
      save
      end
