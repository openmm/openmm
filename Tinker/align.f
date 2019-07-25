c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module align  --  information for structure superposition  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nfit    number of atoms to use in superimposing two structures
c     ifit    atom numbers of pairs of atoms to be superimposed
c     wfit    weights assigned to atom pairs during superposition
c
c
      module align
      implicit none
      integer nfit
      integer, allocatable :: ifit(:,:)
      real*8, allocatable :: wfit(:)
      save
      end
