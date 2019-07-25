c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine version  --  create version number for file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "version" checks the name of a file about to be opened; if
c     if "old" status is passed, the name of the highest current
c     version is returned; if "new" status is passed the filename
c     of the next available unused version is generated
c
c
      subroutine version (filename,status)
      use iounit
      use output
      implicit none
      integer i,leng,trimtext
      integer thousand,hundred
      integer tens,ones
      logical exist
      character*1 digit(0:9)
      character*3 status
      character*240 filename
      character*240 oldfile
      character*240 newfile
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     process the filename and status variables
c
      call lowcase (status)
      leng = trimtext (filename)
c
c     no change is needed if the file doesn't exist
c
      exist = .false.
      if (leng .ne. 0)  inquire (file=filename(1:leng),exist=exist)
      if (.not. exist)  return
c
c     set initial values for the current and next versions
c
      newfile = filename
      oldfile = filename
c
c     append an artificial version number to the filename;
c     currently handles up to 10000 versions of a file
c
      if (.not. noversion) then
         i = 1
         do while (exist)
            i = i + 1
            oldfile = newfile
            thousand = i / 1000
            hundred = (i - 1000*thousand) / 100
            tens = (i - 1000*thousand - 100*hundred) / 10
            ones = i - 1000*thousand - 100*hundred - 10*tens
            if (thousand .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(thousand)
     &                      //digit(hundred)//digit(tens)//digit(ones)
            else if (hundred .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(hundred)
     &                      //digit(tens)//digit(ones)
            else if (tens .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(tens)//digit(ones)
            else
               newfile = filename(1:leng)//'_'//digit(ones)
            end if
            inquire (file=newfile,exist=exist)
         end do
      end if
c
c     set the file name based on the requested status
c
      if (status .eq. 'old') then
         filename = oldfile
      else if (status .eq. 'new') then
         filename = newfile
         inquire (file=filename,exist=exist)
         if (exist) then
            call nextarg (filename,exist)
            if (exist) then
               inquire (file=filename,exist=exist)
            else
               exist = .true.
            end if
            do while (exist)
               write (iout,10)
   10          format (/,' Enter File Name for Coordinate Output :  ',$)
               read (input,20)  filename
   20          format (a240)
               inquire (file=filename,exist=exist)
            end do
         end if
      end if
      return
      end
