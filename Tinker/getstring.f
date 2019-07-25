c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getstring  --  extract double quoted string  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getstring" searches for a quoted text string within an input
c     character string; the region between the first and second
c     double quote is returned as the "text"; if the actual text is
c     too long, only the first part is returned
c
c     variables and parameters:
c
c     string    input character string to be searched
c     text      the quoted text found in the input string
c     next      input with first position of search string;
c                 output with the position following text
c
c
      subroutine getstring (string,text,next)
      use ascii
      implicit none
      integer i,j,k,m
      integer len,length
      integer size,next
      integer code,extent
      integer initial,final
      integer first,last
      integer maxascii
      character*(*) string
      character*(*) text
c
c
c     get the length of input string and output text
c
      length = len(string(next:))
      size = len(text)
c
c     convert first two non-ascii regions to delimiting quotes
c
      maxascii = 126
      initial = next
      final = next + length - 1
      do i = initial, final
         code = ichar(string(i:i))
         if (code .gt. maxascii) then
            string(i:i) = ' '
            do j = i+1, final
               code = ichar(string(j:j))
               if (code .le. maxascii) then
                  string(j-1:j-1) = '"'
                  do k = j+1, final
                     code = ichar(string(k:k))
                     if (code .gt. maxascii) then
                        string(k:k) = '"'
                        do m = k+1, final
                           code = ichar(string(m:m))
                           if (code .gt. maxascii) then
                              string(m:m) = ' '
                           end if
                        end do
                        goto 10
                     end if
                  end do
               end if
            end do
         end if
      end do
   10 continue
c
c     search the string for quoted region of text characters
c
      first = next
      last = 0
      do i = initial, final
         code = ichar(string(i:i))
         if (code .eq. quote) then
            first = i + 1
            do j = first, final
               code = ichar(string(j:j))
               if (code .eq. quote) then
                  last = j - 1
                  next = j + 1
                  goto 20
               end if
            end do
         end if
      end do
   20 continue
c
c     trim the actual word if it is too long to return
c
      extent = last - first + 1
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
      do i = last+1, final
         j = j + 1
         text(j:j) = ' '
      end do
      return
      end
