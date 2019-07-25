c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module domega  --  derivative components over torsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tesum   total energy derivatives over torsions
c     teb     bond stretch derivatives over torsions
c     tea     angle bend derivatives over torsions
c     teba    stretch-bend derivatives over torsions
c     teub    Urey-Bradley derivatives over torsions
c     teaa    angle-angle derivatives over torsions
c     teopb   out-of-plane bend derivatives over torsions
c     teopd   out-of-plane distance derivatives over torsions
c     teid    improper dihedral derivatives over torsions
c     teit    improper torsion derivatives over torsions
c     tet     torsional derivatives over torsions
c     tept    pi-system torsion derivatives over torsions
c     tebt    stretch-torsion derivatives over torsions
c     teat    angle-torsion derivatives over torsions
c     tett    torsion-torsion derivatives over torsions
c     tev     van der Waals derivatives over torsions
c     tec     charge-charge derivatives over torsions
c     tecd    charge-dipole derivatives over torsions
c     ted     dipole-dipole derivatives over torsions
c     tem     atomic multipole derivatives over torsions
c     tep     polarization derivatives over torsions
c     ter     reaction field derivatives over torsions
c     tes     solvation derivatives over torsions
c     telf    metal ligand field derivatives over torsions
c     teg     geometric restraint derivatives over torsions
c     tex     extra energy term derivatives over torsions
c
c
      module domega
      implicit none
      real*8, allocatable :: tesum(:)
      real*8, allocatable :: teb(:)
      real*8, allocatable :: tea(:)
      real*8, allocatable :: teba(:)
      real*8, allocatable :: teub(:)
      real*8, allocatable :: teaa(:)
      real*8, allocatable :: teopb(:)
      real*8, allocatable :: teopd(:)
      real*8, allocatable :: teid(:)
      real*8, allocatable :: teit(:)
      real*8, allocatable :: tet(:)
      real*8, allocatable :: tept(:)
      real*8, allocatable :: tebt(:)
      real*8, allocatable :: teat(:)
      real*8, allocatable :: tett(:)
      real*8, allocatable :: tev(:)
      real*8, allocatable :: tec(:)
      real*8, allocatable :: tecd(:)
      real*8, allocatable :: ted(:)
      real*8, allocatable :: tem(:)
      real*8, allocatable :: tep(:)
      real*8, allocatable :: ter(:)
      real*8, allocatable :: tes(:)
      real*8, allocatable :: telf(:)
      real*8, allocatable :: teg(:)
      real*8, allocatable :: tex(:)
      save
      end
