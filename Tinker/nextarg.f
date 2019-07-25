c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nextarg  --  find next command line argument  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nextarg" finds the next unused command line argument
c     and returns it in the input character string
c
c
      subroutine nextarg (string,exist)
      use argue
      implicit none
      integer i,length
      logical exist
      character*(*) string
c
c
c     initialize the command argument as a blank string
c
      string = '          '
      exist = .false.
c
c     get the next command line argument and mark it as used
c
      if (narg .ne. 0) then
         length = min(len(string),len(arg(maxarg)))
         do i = 1, narg
            if (listarg(i)) then
               listarg(i) = .false.
               string = arg(i)(1:length)
               exist = .true.
               goto 10
            end if
         end do
   10    continue
      end if
      return
      end
