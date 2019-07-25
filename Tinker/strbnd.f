c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module strbnd  --  stretch-bends in current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nstrbnd   total number of stretch-bend interactions
c     isb       angle and bond numbers used in stretch-bend
c     sbk       force constants for stretch-bend terms
c
c
      module strbnd
      implicit none
      integer nstrbnd
      integer, allocatable :: isb(:,:)
      real*8, allocatable :: sbk(:,:)
      save
      end
