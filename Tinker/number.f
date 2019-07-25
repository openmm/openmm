c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  function number  --  convert text string to number  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "number" converts a text numeral into an integer value;
c     the input string must contain only numeric characters
c
c
      function number (string)
      use inform
      use iounit
      implicit none
      integer i,j,number
      integer first,last,trimtext
      integer digit,place(10)
      character*1 letter
      character*(*) string
      data place  / 1, 10, 100, 1000, 10000, 100000, 1000000,
     &              10000000, 100000000, 1000000000 /
c
c
c     initialize the integer value of number to zero
c
      number = 0
c
c     get the first and last nonblank characters
c
      last = trimtext (string)
      if (last .gt. 10) then
         write (iout,10)
   10    format (' NUMBER  --  Input Text String is Too Long')
         return
      end if
      first = 1
      do i = 1, last
         letter = string(i:i)
         if (letter .ne. ' ') then
            first = i
            goto 20
         end if
      end do
   20 continue
c
c     convert the text numeral into an integer number
c
      j = 0
      do i = last, first, -1
         j = j + 1
         letter = string(i:i)
         if (letter .eq. '0') then
            digit = 0
         else if (letter .eq. '1') then
            digit = 1
         else if (letter .eq. '2') then
            digit = 2
         else if (letter .eq. '3') then
            digit = 3
         else if (letter .eq. '4') then
            digit = 4
         else if (letter .eq. '5') then
            digit = 5
         else if (letter .eq. '6') then
            digit = 6
         else if (letter .eq. '7') then
            digit = 7
         else if (letter .eq. '8') then
            digit = 8
         else if (letter .eq. '9') then
            digit = 9
         else
            if (debug) then
               write (iout,30)
   30          format (/,' NUMBER  --  Non-Numeric Characters Found',
     &                    ' in Numeral String')
            end if
            number = 0
            goto 40
         end if
         number = number + digit * place(j)
      end do
   40 continue
      return
      end
