c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine numeral  --  convert number to text string  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "numeral" converts an input integer number into the
c     corresponding right- or left-justified text numeral
c
c     number  integer value of the number to be transformed
c     string  text string to be filled with corresponding numeral
c     size    on input, the minimal acceptable numeral length, if
c               zero then output will be right justified, if
c               nonzero then numeral is left-justified and padded
c               with leading zeros as necessary; upon output, the
c               number of non-blank characters in the numeral
c
c
      subroutine numeral (number,string,size)
      implicit none
      integer i,number,size
      integer multi,pos
      integer length,minsize,len
      integer million,hunthou
      integer tenthou,thousand
      integer hundred,tens,ones
      logical right,negative
      character*1 digit(0:9)
      character*(*) string
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     set justification and size bounds for numeral string
c
      if (size .eq. 0) then
         right = .true.
         size = 1
      else
         right = .false.
      end if
      minsize = size
      length = len(string)
c
c     test the sign of the original number
c
      if (number .ge. 0) then
         negative = .false.
      else
         negative = .true.
         number = -number
      end if
c
c     use modulo arithmetic to find place-holding digits
c
      million = number / 1000000
      multi = 1000000 * million
      hunthou = (number-multi) / 100000
      multi = multi + 100000*hunthou
      tenthou = (number-multi) / 10000
      multi = multi + 10000*tenthou
      thousand = (number-multi) / 1000
      multi = multi + 1000*thousand
      hundred = (number-multi) / 100
      multi = multi + 100*hundred
      tens = (number-multi) / 10
      multi = multi + 10*tens
      ones = number - multi
c
c     find the correct length to be used for the numeral
c
      if (million .ne. 0) then
         size = 7
      else if (hunthou .ne. 0) then
         size = 6
      else if (tenthou .ne. 0) then
         size = 5
      else if (thousand .ne. 0) then
         size = 4
      else if (hundred .ne. 0) then
         size = 3
      else if (tens .ne. 0) then
         size = 2
      else
         size = 1
      end if
      size = min(size,length)
      size = max(size,minsize)
c
c     convert individual digits to a string of numerals
c
      if (size .eq. 7) then
         string(1:1) = digit(million)
         string(2:2) = digit(hunthou)
         string(3:3) = digit(tenthou)
         string(4:4) = digit(thousand)
         string(5:5) = digit(hundred)
         string(6:6) = digit(tens)
         string(7:7) = digit(ones)
      else if (size .eq. 6) then
         string(1:1) = digit(hunthou)
         string(2:2) = digit(tenthou)
         string(3:3) = digit(thousand)
         string(4:4) = digit(hundred)
         string(5:5) = digit(tens)
         string(6:6) = digit(ones)
      else if (size .eq. 5) then
         string(1:1) = digit(tenthou)
         string(2:2) = digit(thousand)
         string(3:3) = digit(hundred)
         string(4:4) = digit(tens)
         string(5:5) = digit(ones)
      else if (size .eq. 4) then
         string(1:1) = digit(thousand)
         string(2:2) = digit(hundred)
         string(3:3) = digit(tens)
         string(4:4) = digit(ones)
      else if (size .eq. 3) then
         string(1:1) = digit(hundred)
         string(2:2) = digit(tens)
         string(3:3) = digit(ones)
      else if (size .eq. 2) then
         string(1:1) = digit(tens)
         string(2:2) = digit(ones)
      else
         string(1:1) = digit(ones)
      end if
c
c     right-justify if desired, and pad with blanks
c
      if (right) then
         do i = size, 1, -1
            pos = length - size + i
            string(pos:pos) = string(i:i)
         end do
         do i = 1, length-size
            string(i:i) = ' '
         end do
      else
         do i = size+1, length
            string(i:i) = ' '
         end do
      end if
      return
      end
