c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module precis  --  machine precision tolerance values  ##
c     ##                                                         ##
c     #############################################################
c
c
c     tiny    the smallest positive floating point value
c     small   the smallest relative floating point spacing
c     huge    the largest relative floating point spacing
c
c
      module precis
      implicit none
      real*8 tiny
      real*8 small
      real*8 huge
      save
      end
