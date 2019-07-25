c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module piorbs  --  conjugated system in current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     norbit    total number of pisystem orbitals in the system
c     nconj     total number of separate conjugated piystems
c     reorbit   number of evaluations between orbital updates
c     nbpi      total number of bonds affected by the pisystem
c     ntpi      total number of torsions affected by the pisystem
c     iorbit    numbers of the atoms containing pisystem orbitals
c     iconj     first and last atom of each pisystem in the list
c     kconj     contiguous list of atoms in each pisystem
c     piperp    atoms defining a normal plane to each orbital
c     ibpi      bond and piatom numbers for each pisystem bond
c     itpi      torsion and pibond numbers for each pisystem torsion
c     pbpl      pi-bond orders for bonds in "planar" pisystem
c     pnpl      pi-bond orders for bonds in "nonplanar" pisystem
c     listpi    atom list indicating whether each atom has an orbital
c
c
      module piorbs
      implicit none
      integer norbit
      integer nconj
      integer reorbit
      integer nbpi
      integer ntpi
      integer, allocatable :: iorbit(:)
      integer, allocatable :: iconj(:,:)
      integer, allocatable :: kconj(:)
      integer, allocatable :: piperp(:,:)
      integer, allocatable :: ibpi(:,:)
      integer, allocatable :: itpi(:,:)
      real*8, allocatable :: pbpl(:)
      real*8, allocatable :: pnpl(:)
      logical, allocatable :: listpi(:)
      save
      end
