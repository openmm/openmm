c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module syntrn  --  synchronous transit path definition  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tpath   value of the path coordinate (0=reactant, 1=product)
c     ppath   path coordinate for extra point in quadratic transit
c     xmin1   reactant coordinates as array of optimization variables
c     xmin2   product coordinates as array of optimization variables
c     xm      extra coordinate set for quadratic synchronous transit
c
c
      module syntrn
      implicit none
      real*8 tpath
      real*8 ppath
      real*8, allocatable :: xmin1(:)
      real*8, allocatable :: xmin2(:)
      real*8, allocatable :: xm(:)
      save
      end
