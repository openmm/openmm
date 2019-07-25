c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine basefile  --  get base prefix from a filename  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "basefile" extracts from an input filename the portion
c     consisting of any directory name and the base filename;
c     also reads any keyfile and sets information level values
c
c
      subroutine basefile (string)
      use ascii
      use files
      implicit none
      integer i,k,trimtext
      character*1 letter
      character*240 string
      character*240 prefix
c
c
c     account for home directory abbreviation in filename
c
      if (string(1:2) .eq. '~/') then
         call getenv ('HOME',prefix)
         string = prefix(1:trimtext(prefix))//
     &               string(2:trimtext(string))
      end if
c
c     store the input filename and find its full length
c
      filename = string
      leng = trimtext (string)
c
c     count the number of characters prior to any extension
c
      k = leng
      do i = 1, leng
         letter = string(i:i)
         if (letter .eq. '/')  k = leng
c        if (letter .eq. '\')  k = leng
         if (ichar(letter) .eq. backslash)  k = leng
         if (letter .eq. ']')  k = leng
         if (letter .eq. ':')  k = leng
         if (letter .eq. '~')  k = leng
         if (letter .eq. '.')  k = i - 1
      end do
      leng = min(leng,k)
c
c     find the length of any directory name prefix
c
      k = 0
      do i = leng, 1, -1
         letter = string(i:i)
         if (letter .eq. '/')  k = i
c        if (letter .eq. '\')  k = i
         if (ichar(letter) .eq. backslash)  k = i
         if (letter .eq. ']')  k = i
         if (letter .eq. ':')  k = i
         if (letter .eq. '~')  k = i
c        if (letter .eq. '.')  k = i
         if (k .ne. 0)  goto 10
      end do
   10 continue
      ldir = k
c
c     read and store the keywords from the keyfile
c
      call getkey
c
c     get the information level and output style
c
      call control
      return
      end
