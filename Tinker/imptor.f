c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module imptor  --  improper torsions in current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nitors   total number of improper torsional angles in the system
c     iitors   numbers of the atoms in each improper torsional angle
c     itors1   1-fold amplitude and phase for each improper torsion
c     itors2   2-fold amplitude and phase for each improper torsion
c     itors3   3-fold amplitude and phase for each improper torsion
c
c
      module imptor
      implicit none
      integer nitors
      integer, allocatable :: iitors(:,:)
      real*8, allocatable :: itors1(:,:)
      real*8, allocatable :: itors2(:,:)
      real*8, allocatable :: itors3(:,:)
      save
      end
