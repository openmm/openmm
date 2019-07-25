c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module rgddyn  --  rigid body MD velocities and momenta  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xcmo    x-component from each atom to center of rigid body
c     ycmo    y-component from each atom to center of rigid body
c     zcmo    z-component from each atom to center of rigid body
c     vcm     current translational velocity of each rigid body
c     wcm     current angular velocity of each rigid body
c     lm      current angular momentum of each rigid body
c     vc      half-step translational velocity for kinetic energy
c     wc      half-step angular velocity for kinetic energy
c     linear  logical flag to mark group as linear or nonlinear
c
c
      module rgddyn
      implicit none
      real*8, allocatable :: xcmo(:)
      real*8, allocatable :: ycmo(:)
      real*8, allocatable :: zcmo(:)
      real*8, allocatable :: vcm(:,:)
      real*8, allocatable :: wcm(:,:)
      real*8, allocatable :: lm(:,:)
      real*8, allocatable :: vc(:,:)
      real*8, allocatable :: wc(:,:)
      logical, allocatable :: linear(:)
      save
      end
