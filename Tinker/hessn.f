c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module hessn  --  Cartesian Hessian elements for one atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     hessx   Hessian elements for x-component of current atom
c     hessy   Hessian elements for y-component of current atom
c     hessz   Hessian elements for z-component of current atom
c
c
      module hessn
      implicit none
      real*8, allocatable :: hessx(:,:)
      real*8, allocatable :: hessy(:,:)
      real*8, allocatable :: hessz(:,:)
      save
      end
