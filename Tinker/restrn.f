c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module restrn  --  parameters for geometrical restraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     npfix      number of position restraints to be applied
c     ndfix      number of distance restraints to be applied
c     nafix      number of angle restraints to be applied
c     ntfix      number of torsional restraints to be applied
c     ngfix      number of group distance restraints to be applied
c     nchir      number of chirality restraints to be applied
c     ipfix      atom number involved in each position restraint
c     kpfix      flags to use x-, y-, z-coordinate position restraints
c     idfix      atom numbers defining each distance restraint
c     iafix      atom numbers defining each angle restraint
c     itfix      atom numbers defining each torsional restraint
c     igfix      group numbers defining each group distance restraint
c     ichir      atom numbers defining each chirality restraint
c     depth      depth of shallow Gaussian basin restraint
c     width      exponential width coefficient of Gaussian basin
c     rwall      radius of spherical droplet boundary restraint
c     xpfix      x-coordinate target for each restrained position
c     ypfix      y-coordinate target for each restrained position
c     zpfix      z-coordinate target for each restrained position
c     pfix       force constant and flat-well range for each position
c     dfix       force constant and target range for each distance
c     afix       force constant and target range for each angle
c     tfix       force constant and target range for each torsion
c     gfix       force constant and target range for each group distance
c     chir       force constant and target range for chiral centers
c     use_basin  logical flag governing use of Gaussian basin
c     use_wall   logical flag governing use of droplet boundary
c
c
      module restrn
      implicit none
      integer npfix
      integer ndfix
      integer nafix
      integer ntfix
      integer ngfix
      integer nchir
      integer, allocatable :: ipfix(:)
      integer, allocatable :: kpfix(:,:)
      integer, allocatable :: idfix(:,:)
      integer, allocatable :: iafix(:,:)
      integer, allocatable :: itfix(:,:)
      integer, allocatable :: igfix(:,:)
      integer, allocatable :: ichir(:,:)
      real*8 depth
      real*8 width
      real*8 rwall
      real*8, allocatable :: xpfix(:)
      real*8, allocatable :: ypfix(:)
      real*8, allocatable :: zpfix(:)
      real*8, allocatable :: pfix(:,:)
      real*8, allocatable :: dfix(:,:)
      real*8, allocatable :: afix(:,:)
      real*8, allocatable :: tfix(:,:)
      real*8, allocatable :: gfix(:,:)
      real*8, allocatable :: chir(:,:)
      logical use_basin
      logical use_wall
      save
      end
