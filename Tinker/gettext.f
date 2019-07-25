c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine gettext  --  extract text from a string  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "gettext" searches an input string for the first string of
c     non-blank characters; the region from a non-blank character
c     to the first space or tab is returned as "text"; if the
c     actual text is too long, only the first part is returned
c
c     variables and parameters:
c
c     string    input character string to be searched
c     text      output with the first text string found
c     next      input with first position of search string;
c                 output with the position following text
c
c
      subroutine gettext (string,text,next)
      use ascii
      implicit none
      integer i,j
      integer len,length
      integer size,next
      integer first,last
      integer code,extent
      integer initial,final
      character*(*) string
      character*(*) text
c
c
c     get the length of input string and output text
c
      length = len(string(next:))
      size = len(text)
c
c     search the string for the first non-blank character
c
      first = next
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         code = ichar(string(i:i))
         if (code.ne.space .and. code.ne.tab) then
            first = i
            do j = i+1, final
               code = ichar(string(j:j))
               if (code.eq.space .or. code.eq.tab) then
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
c     trim the actual text if it is too long to return
c
      extent = next - first
      final = first + size - 1
      if (extent .gt. size)  last = final
c
c     transfer the text into the return string
c
      j = 0
      do i = first, last
         j = j + 1
         text(j:j) = string(i:i)
      end do
      do i = next, final
         j = j + 1
         text(j:j) = ' '
      end do
      return
      end
