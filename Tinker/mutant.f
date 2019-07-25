c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module mutant  --  free energy calculation hybrid atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for van der Waals potentials
c     elambda    state weighting value for electrostatic potentials
c     scexp      scale factor for soft core buffered 14-7 potential
c     scalpha    scale factor for soft core buffered 14-7 potential
c     mut        true if an atom is to be mutated, false otherwise
c
c
      module mutant
      implicit none
      integer nmut
      integer, allocatable :: imut(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: class0(:)
      integer, allocatable :: type1(:)
      integer, allocatable :: class1(:)
      real*8 lambda
      real*8 vlambda
      real*8 elambda
      real*8 scexp
      real*8 scalpha
      logical, allocatable :: mut(:)
      save
      end
