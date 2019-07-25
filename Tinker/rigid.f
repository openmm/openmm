c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module rigid  --  rigid body coordinates for atom groups  ##
c     ##                                                            ##
c     ################################################################
c
c
c     xrb         rigid body reference x-coordinate for each atom
c     yrb         rigid body reference y-coordinate for each atom
c     zrb         rigid body reference z-coordinate for each atom
c     rbc         current rigid body coordinates for each group
c     use_rigid   flag to mark use of rigid body coordinate system
c
c
      module rigid
      implicit none
      real*8, allocatable :: xrb(:)
      real*8, allocatable :: yrb(:)
      real*8, allocatable :: zrb(:)
      real*8, allocatable :: rbc(:,:)
      logical use_rigid
      save
      end
