c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module orbits  --  conjugated pisystem orbital energies  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     qorb    number of pi-electrons contributed by each atom
c     worb    ionization potential of each pisystem atom
c     emorb   repulsion integral for each pisystem atom
c
c
      module orbits
      implicit none
      real*8, allocatable :: qorb(:)
      real*8, allocatable :: worb(:)
      real*8, allocatable :: emorb(:)
      save
      end
