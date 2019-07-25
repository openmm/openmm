c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prterr  --  output coordinates upon error  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prterr" writes out a set of coordinates to a disk
c     file prior to aborting on a serious error
c
c
      subroutine prterr
      use files
      use output
      implicit none
      integer ierr,freeunit
      character*240 errorfile
c
c
c     write the current coordinates to a file after an error
c
      ierr = freeunit ()
      errorfile = filename(1:leng)//'.err'
      call version (errorfile,'new')
      open (unit=ierr,file=errorfile,status='new')
      if (coordtype .eq. 'CARTESIAN') then
         call prtxyz (ierr)
      else if (coordtype .eq. 'INTERNAL') then
         call prtint (ierr)
      else if (coordtype .eq. 'RIGIDBODY') then
         call prtxyz (ierr)
      end if
      close (unit=ierr)
      return
      end
