c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module analyz  --  energy components partitioned to atoms  ##
c     ##                                                             ##
c     #################################################################
c
c
c     aesum   total potential energy partitioned over atoms
c     aeb     bond stretch energy partitioned over atoms
c     aea     angle bend energy partitioned over atoms
c     aeba    stretch-bend energy partitioned over atoms
c     aeub    Urey-Bradley energy partitioned over atoms
c     aeaa    angle-angle energy partitioned over atoms
c     aeopb   out-of-plane bend energy partitioned over atoms
c     aeopd   out-of-plane distance energy partitioned over atoms
c     aeid    improper dihedral energy partitioned over atoms
c     aeit    improper torsion energy partitioned over atoms
c     aet     torsional energy partitioned over atoms
c     aept    pi-system torsion energy partitioned over atoms
c     aebt    stretch-torsion energy partitioned over atoms
c     aeat    angle-torsion energy partitioned over atoms
c     aett    torsion-torsion energy partitioned over atoms
c     aev     van der Waals energy partitioned over atoms
c     act     charge transfer energy partitioned over atoms
c     aec     charge-charge energy partitioned over atoms
c     aecd    charge-dipole energy partitioned over atoms
c     aed     dipole-dipole energy partitioned over atoms
c     aem     multipole energy partitioned over atoms
c     aep     polarization energy partitioned over atoms
c     aer     reaction field energy partitioned over atoms
c     aes     solvation energy partitioned over atoms
c     aelf    metal ligand field energy partitioned over atoms
c     aeg     geometric restraint energy partitioned over atoms
c     aex     extra energy term partitioned over atoms
c
c
      module analyz
      implicit none
      real*8, allocatable :: aesum(:)
      real*8, allocatable :: aeb(:)
      real*8, allocatable :: aea(:)
      real*8, allocatable :: aeba(:)
      real*8, allocatable :: aeub(:)
      real*8, allocatable :: aeaa(:)
      real*8, allocatable :: aeopb(:)
      real*8, allocatable :: aeopd(:)
      real*8, allocatable :: aeid(:)
      real*8, allocatable :: aeit(:)
      real*8, allocatable :: aet(:)
      real*8, allocatable :: aept(:)
      real*8, allocatable :: aebt(:)
      real*8, allocatable :: aeat(:)
      real*8, allocatable :: aett(:)
      real*8, allocatable :: aev(:)
CW!!!
      real*8, allocatable :: aect(:)
      real*8, allocatable :: aec(:)
      real*8, allocatable :: aecd(:)
      real*8, allocatable :: aed(:)
      real*8, allocatable :: aem(:)
      real*8, allocatable :: aep(:)
      real*8, allocatable :: aer(:)
      real*8, allocatable :: aes(:)
      real*8, allocatable :: aelf(:)
      real*8, allocatable :: aeg(:)
      real*8, allocatable :: aex(:)
      save
      end
