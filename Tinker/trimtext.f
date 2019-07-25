c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function trimtext  --  find last non-blank character  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "trimtext" finds and returns the location of the last
c     non-blank character before the first null character in
c     an input text string; the function returns zero if no
c     such character is found
c
c
      function trimtext (string)
      implicit none
      integer i,size,last
      integer len,trimtext
      character*1 char
      character*1 null
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first null character
c
      trimtext = 0
      size = len(string)
      null = char(0)
      last = size
      do i = 1, size
         if (string(i:i) .eq. null) then
            last = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      do i = last, 1, -1
         if (string(i:i) .gt. ' ') then
            trimtext = i
            goto 20
         end if
      end do
   20 continue
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine trimhead  --  remove spaces before first text  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "trimhead" removes blank spaces before the first non-blank
c     character in a text string by shifting the string to the left
c
c
      subroutine trimhead (string)
      implicit none
      integer i,j,k
      character*240 string
      character*240 temp
c
c
c     loop over characters, removing blank beginning spaces
c
      do i = 1, 240
         temp(i:i) = ' '
      end do
      j = 0
      k = 0
      do i = 1, 240
         if (string(i:i) .ne. ' ')  j = 1
         if (j .eq. 1) then
            k = k + 1
            temp(k:k) = string(i:i)
         end if
      end do
      do i = 1, 240
         string(i:i) = temp(i:i)
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine justify  --  convert string to right justified  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "justify" converts a text string to right justified format
c     with leading blank spaces
c
c
      subroutine justify (string)
      implicit none
      integer i,k,len
      integer size,last
      character*1 char
      character*1 null
      character*1 letter
      character*(*) string
c
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      size = len(string)
      null = char(0)
      last = 0
      do i = size, 1, -1
         letter = string(i:i)
         if (letter.ne.' ' .and. letter.ne.null) then
            last = i
            goto 10
         end if
      end do
   10 continue
c
c     move string to the right and pad with leading blanks
c
      do i = last, 1, -1
         k = i + size - last
         string(k:k) = string(i:i)
      end do
      do i = 1, size-last
         string(i:i) = ' '
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine upcase  --  convert string to all upper case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "upcase" converts a text string to all upper case letters
c
c
      subroutine upcase (string)
      implicit none
      integer i,size,len
      integer code,ichar
      character*1 char
      character*1 letter
      character*(*) string
c
c
c     convert lower case to upper case one letter at a time
c
      size = len(string)
      do i = 1, size
         letter = string(i:i)
         code = ichar(letter)
         if (letter.ge.'a' .and. letter.le.'z')
     &      string(i:i) = char(code-32)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lowcase  --  convert string to all lower case  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lowcase" converts a text string to all lower case letters
c
c
      subroutine lowcase (string)
      implicit none
      integer i,size
      integer code,ichar
      character*1 char
      character*1 letter
      character*(*) string
c
c
c     convert upper case to lower case one letter at a time
c
      size = len(string)
      do i = 1, size
         letter = string(i:i)
         code = ichar(letter)
         if (letter.ge.'A' .and. letter.le.'Z')
     &      string(i:i) = char(code+32)
      end do
      return
      end
