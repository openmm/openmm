c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module angbnd  --  bond angle bends in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nangle   total number of angle bends in the system
c     iang     numbers of the atoms in each angle bend
c     ak       harmonic angle force constant (kcal/mole/rad**2)
c     anat     ideal bond angle or phase shift angle (degrees)
c     afld     periodicity for Fourier angle bending term
c
c
      module angbnd
      implicit none
      integer nangle
      integer, allocatable :: iang(:,:)
      real*8, allocatable :: ak(:)
      real*8, allocatable :: anat(:)
      real*8, allocatable :: afld(:)
      save
      end
