c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine fatal  --  terminate the program abnormally  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "fatal" terminates execution due to a user request, a severe
c     error or some other nonstandard condition
c
c
      subroutine fatal
      use iounit
      implicit none
c
c
c     print a final warning message, then do final cleanup
c
      write (iout,10)
   10 format (/,' TINKER is Unable to Continue; Terminating',
     &           ' the Current Calculation',/)
      call final
c
c     exit is not standard Fortran, but stop can give errors
c
      call exit
c     stop
      end
