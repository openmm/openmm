c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine suffix  --  test for default file extension  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "suffix" checks a filename for the presence of an extension,
c     and appends an extension and version if none is found
c
c
      subroutine suffix (filename,extension,status)
      use ascii
      implicit none
      integer i,leng,lext
      integer last,trimtext
      logical exist
      character*1 letter
      character*3 status
      character*(*) filename
      character*(*) extension
c
c
c     get the length of the current filename
c
      leng = trimtext (filename)
      lext = trimtext (extension)
c
c     check for an extension on the current filename
c
      last = leng
      do i = 1, leng
         letter = filename(i:i)
         if (letter .eq. '/')  last = leng
c        if (letter .eq. '\')  last = leng
         if (ichar(letter) .eq. backslash)  last = leng
         if (letter .eq. ']')  last = leng
         if (letter .eq. ':')  last = leng
         if (letter .eq. '~')  last = leng
         if (letter .eq. '.')  last = i - 1
      end do
c
c     append an extension or version as appropriate
c
      if (last .eq. leng) then
         exist = .false.
         if (leng .ne. 0) then
            inquire (file=filename(1:leng),exist=exist)
         end if
         if (.not. exist) then
            filename = filename(1:leng)//'.'//extension(1:lext)
            call version (filename,status)
         end if
      else if (status .eq. 'new') then
         call version (filename,status)
      end if
      return
      end
