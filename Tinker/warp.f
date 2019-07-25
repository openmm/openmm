c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module warp  --  potential surface smoothing parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     deform       value of the smoothing deformation parameter
c     difft        diffusion coefficient for torsional potential
c     diffv        diffusion coefficient for van der Waals potential
c     diffc        diffusion coefficient for charge-charge potential
c     m2           second moment of the GDA gaussian for each atom
c     use_smooth   flag to use a potential energy smoothing method
c     use_dem      flag to use diffusion equation method potential
c     use_gda      flag to use gaussian density annealing potential
c     use_tophat   flag to use analytical tophat smoothed potential
c     use_stophat  flag to use shifted tophat smoothed potential
c
c
      module warp
      implicit none
      real*8 deform
      real*8 difft
      real*8 diffv
      real*8 diffc
      real*8, allocatable :: m2(:)
      logical use_smooth
      logical use_dem
      logical use_gda
      logical use_tophat
      logical use_stophat
      save
      end
