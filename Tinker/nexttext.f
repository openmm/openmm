c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function nexttext  --  find next non-blank character  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "nexttext" finds and returns the location of the first
c     non-blank character within an input text string; zero
c     is returned if no such character is found
c
c
      function nexttext (string)
      implicit none
      integer i,size
      integer len,nexttext
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first non-blank character
c
      nexttext = 0
      size = len(string)
      do i = 1, size
         if (string(i:i) .gt. ' ') then
            nexttext = i
            goto 10
         end if
      end do
   10 continue
      return
      end
