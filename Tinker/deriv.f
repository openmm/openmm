c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module deriv  --  Cartesian coord derivative components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-system torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     deat    angle-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     dect    charge transfer Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decd    charge-dipole Cartesian coordinate derivatives
c     ded     dipole-dipole Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     der     reaction field Cartesian coordinate derivatives
c     des     solvation Cartesian coordinate derivatives
c     delf    metal ligand field Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c
c
      module deriv
      implicit none
      real*8, allocatable :: desum(:,:)
      real*8, allocatable :: deb(:,:)
      real*8, allocatable :: dea(:,:)
      real*8, allocatable :: deba(:,:)
      real*8, allocatable :: deub(:,:)
      real*8, allocatable :: deaa(:,:)
      real*8, allocatable :: deopb(:,:)
      real*8, allocatable :: deopd(:,:)
      real*8, allocatable :: deid(:,:)
      real*8, allocatable :: deit(:,:)
      real*8, allocatable :: det(:,:)
      real*8, allocatable :: dept(:,:)
      real*8, allocatable :: debt(:,:)
      real*8, allocatable :: deat(:,:)
      real*8, allocatable :: dett(:,:)
      real*8, allocatable :: dev(:,:)
      real*8, allocatable :: dect(:,:)
      real*8, allocatable :: dec(:,:)
      real*8, allocatable :: decd(:,:)
      real*8, allocatable :: ded(:,:)
      real*8, allocatable :: dem(:,:)
      real*8, allocatable :: dep(:,:)
      real*8, allocatable :: der(:,:)
      real*8, allocatable :: des(:,:)
      real*8, allocatable :: delf(:,:)
      real*8, allocatable :: deg(:,:)
      real*8, allocatable :: dex(:,:)
      save
      end
