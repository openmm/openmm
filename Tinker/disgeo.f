c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module disgeo  --  distance geometry bounds & parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     vdwmax      maximum value of hard sphere sum for an atom pair
c     compact     index of local distance compaction on embedding
c     pathmax     maximum value of upper bound after smoothing
c     dbnd        distance geometry upper and lower bounds matrix
c     georad      hard sphere radii for distance geometry atoms
c     use_invert  flag to use enantiomer closest to input structure
c     use_anneal  flag to use simulated annealing refinement
c
c
      module disgeo
      implicit none
      real*8 vdwmax
      real*8 compact
      real*8 pathmax
      real*8, allocatable :: dbnd(:,:)
      real*8, allocatable :: georad(:)
      logical use_invert
      logical use_anneal
      save
      end
