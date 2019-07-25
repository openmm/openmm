c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module opbend  --  out-of-plane bends in current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     nopbend   total number of out-of-plane bends in the system
c     iopb      bond angle numbers used in out-of-plane bending
c     opbk      force constant values for out-of-plane bending
c
c
      module opbend
      implicit none
      integer nopbend
      integer, allocatable :: iopb(:)
      real*8, allocatable :: opbk(:)
      save
      end
