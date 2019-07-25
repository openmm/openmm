c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module dipole  --  bond dipoles in current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ndipole   total number of dipoles in the system
c     idpl      numbers of atoms that define each dipole
c     bdpl      magnitude of each of the dipoles (Debyes)
c     sdpl      position of each dipole between defining atoms
c
c
      module dipole
      implicit none
      integer ndipole
      integer, allocatable :: idpl(:,:)
      real*8, allocatable :: bdpl(:)
      real*8, allocatable :: sdpl(:)
      save
      end
