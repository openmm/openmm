c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getword  --  extract first word from a string  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getword" searches an input string for the first alphabetic
c     character (A-Z or a-z); the region from this first character
c     to the first blank space or separator is returned as a "word";
c     if the actual word is too long, only the first part is returned
c
c     variables and parameters:
c
c     string    input character string to be searched
c     word      output with the first word in the string
c     next      input with first position of search string;
c                 output with the position following word
c
c
      subroutine getword (string,word,next)
      use ascii
      implicit none
      integer i,j
      integer len,length
      integer size,next
      integer first,last
      integer code,extent
      integer initial,final
      character*1 letter
      character*(*) string
      character*(*) word
c
c
c     get the length of input string and output word
c
      length = len(string(next:))
      size = len(word)
c
c     search the string for the first alphabetic character
c
      first = next
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         letter = string(i:i)
         if ((letter.ge.'A' .and. letter.le.'Z') .or.
     &       (letter.ge.'a' .and. letter.le.'z')) then
            first = i
            do j = i+1, final
               code = ichar(string(j:j))
               if (code.eq.space .or. code.eq.tab .or.
     &             code.eq.comma .or. code.eq.colon .or.
     &             code.eq.semicolon) then
                  last = j - 1
                  next = j
                  goto 10
               end if
            end do
            last = final
            next = last + 1
         end if
      end do
   10 continue
c
c     trim the actual word if it is too long to return
c
      extent = next - first
      final = first + size - 1
      if (extent .gt. size)  last = final
c
c     transfer the word into the return string
c
      j = 0
      do i = first, last
         j = j + 1
         word(j:j) = string(i:i)
      end do
      do i = next, final
         j = j + 1
         word(j:j) = ' '
      end do
c
c     skip over the next character when it is a separator
c
      code = ichar(string(next:next))
      if (code.eq.tab .or. code.eq.comma .or.
     &    code.eq.colon .or. code.eq.semicolon) then
          next = next + 1
      end if
      return
      end
