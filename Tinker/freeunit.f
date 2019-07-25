c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function freeunit  --  gets an unopened logical unit  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "freeunit" finds an unopened Fortran I/O unit and returns
c     its numerical value from 1 to 99; the units already assigned
c     to "input" and "iout" (usually 5 and 6) are skipped since
c     they have special meaning as the default I/O units
c
c
      function freeunit ()
      use iounit
      implicit none
      integer freeunit
      logical used
c
c
c     try each logical unit until an unopened one is found
c
      freeunit = 0
      used = .true.
      do while (used)
         freeunit = freeunit + 1
         if (freeunit.ne.input .and. freeunit.ne.iout) then
            if (freeunit .gt. 99) then
               write (iout,10)
   10          format (/,' FREEUNIT  --  No Available Fortran',
     &                    ' I/O Units')
               call fatal
            end if
            inquire (unit=freeunit,opened=used)
         end if
      end do
      return
      end
