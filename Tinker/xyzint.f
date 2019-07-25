c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program xyzint  --  Cartesian to internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "xyzint" takes as input a Cartesian coordinates file, then
c     converts to and writes out an internal coordinates file
c
c
      program xyzint
      use files
      use iounit
      use titles
      implicit none
      integer izmt,mode
      integer next,freeunit
      logical exist
      character*1 answer
      character*240 intfile
      character*240 record
c
c
c     get and read the Cartesian coordinates file
c
      call initial
      call getxyz
      write (iout,10)  title(1:ltitle)
   10 format (/,' Title :  ',a)
c
c     set the mode for conversion to internal coordinates
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,20)
   20    format (/,' Template (T), Dihedrals (D), Manual (M)',
     &              ' or Automatic [A] :  ',$)
         read (input,30)  record
   30    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      mode = 0
      if (answer .eq. 'M') then
         mode = 1
      else if (answer .eq. 'T') then
         mode = 2
         intfile = filename(1:leng)//'.int'
         call version (intfile,'old')
         inquire (file=intfile,exist=exist)
         if (exist) then
            izmt = freeunit ()
            open (unit=izmt,file=intfile,status='old')
            rewind (unit=izmt)
            call readint (izmt)
            close (unit=izmt)
         else
            mode = 0
            write (iout,40)
   40       format (/,' XYZINT  --  Warning, Template File was',
     &                 ' not Found')
         end if
      else if (answer .eq. 'D') then
         mode = 3
      end if
c
c     convert from Cartesian to internal coordinates
c
      call makeint (mode)
c
c     write out the internal coordinates file
c
      izmt = freeunit ()
      intfile = filename(1:leng)//'.int'
      call version (intfile,'new')
      open (unit=izmt,file=intfile,status='new')
      call prtint (izmt)
      close (unit=izmt)
c
c     perform any final tasks before program exit
c
      call final
      end
