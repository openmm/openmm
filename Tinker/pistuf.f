c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module pistuf  --  bond order-related pisystem parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     bkpi     bond stretch force constants for pi-bond order of 1.0
c     blpi     ideal bond length values for a pi-bond order of 1.0
c     kslope   rate of force constant decrease with bond order decrease
c     lslope   rate of bond length increase with a bond order decrease
c     torsp2   2-fold torsional energy barrier for pi-bond order of 1.0
c
c
      module pistuf
      implicit none
      real*8, allocatable :: bkpi(:)
      real*8, allocatable :: blpi(:)
      real*8, allocatable :: kslope(:)
      real*8, allocatable :: lslope(:)
      real*8, allocatable :: torsp2(:)
      save
      end
