c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module light  --  method of lights pair neighbors indices  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nlight  total number of sites for method of lights calculation
c     kbx     low index of neighbors of each site in the x-sorted list
c     kby     low index of neighbors of each site in the y-sorted list
c     kbz     low index of neighbors of each site in the z-sorted list
c     kex     high index of neighbors of each site in the x-sorted list
c     key     high index of neighbors of each site in the y-sorted list
c     kez     high index of neighbors of each site in the z-sorted list
c     locx    maps the x-sorted list into original interaction list
c     locy    maps the y-sorted list into original interaction list
c     locz    maps the z-sorted list into original interaction list
c     rgx     maps the original interaction list into x-sorted list
c     rgy     maps the original interaction list into y-sorted list
c     rgz     maps the original interaction list into z-sorted list
c
c
      module light
      implicit none
      integer nlight
      integer, allocatable :: kbx(:)
      integer, allocatable :: kby(:)
      integer, allocatable :: kbz(:)
      integer, allocatable :: kex(:)
      integer, allocatable :: key(:)
      integer, allocatable :: kez(:)
      integer, allocatable :: locx(:)
      integer, allocatable :: locy(:)
      integer, allocatable :: locz(:)
      integer, allocatable :: rgx(:)
      integer, allocatable :: rgy(:)
      integer, allocatable :: rgz(:)
      save
      end
