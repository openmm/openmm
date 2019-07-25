c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  function precise  --  determine machine precision values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "precise" finds a machine precision value as selected by the
c     input argument: (1) the smallest positive floating point value,
c     (2) the smallest relative floating point spacing, or (3) the
c     largest relative floating point spacing
c
c
      function precise (i)
      implicit none
      integer i
      real*8 precise,value
      real*8 zero,one,delta
c
c
c     set values for zero, one and multiplicative factor
c
      zero = 0.0d0
      one = 1.0d0
      delta = 1.1d0
      precise = one
c
c     find the smallest positive floating point value;
c     hard coded minimum of 0.24x10-307 is a patch needed
c     to avoid an infinite loop on some machine types
c
      if (i .eq. 1) then
c        do while (precise .ne. zero)
         do while (precise .ge. 0.24d-307)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the smallest relative floating point spacing
c
      else if (i .eq. 2) then
         do while (one+precise .ne. one)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the largest relative floating point spacing
c
      else if (i .eq. 3) then
         do while (one+precise .ne. precise)
           value = precise
           precise = precise * delta
         end do
         precise = value
      end if
      return
      end
