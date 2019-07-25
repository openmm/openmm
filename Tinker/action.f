c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module action  --  total number of each energy term type  ##
c     ##                                                            ##
c     ################################################################
c
c
c     neb     number of bond stretch energy terms computed
c     nea     number of angle bend energy terms computed
c     neba    number of stretch-bend energy terms computed
c     neub    number of Urey-Bradley energy terms computed
c     neaa    number of angle-angle energy terms computed
c     neopb   number of out-of-plane bend energy terms computed
c     neopd   number of out-of-plane distance energy terms computed
c     neid    number of improper dihedral energy terms computed
c     neit    number of improper torsion energy terms computed
c     net     number of torsional energy terms computed
c     nept    number of pi-system torsion energy terms computed
c     nebt    number of stretch-torsion energy terms computed
c     neat    number of angle-torsion energy terms computed
c     nett    number of torsion-torsion energy terms computed
c     nev     number of van der Waals energy terms computed
c     nect     number of van der Waals energy terms computed
c     nec     number of charge-charge energy terms computed
c     necd    number of charge-dipole energy terms computed
c     ned     number of dipole-dipole energy terms computed
c     nem     number of multipole energy terms computed
c     nep     number of polarization energy terms computed
c     new     number of Ewald summation energy terms computed
c     ner     number of reaction field energy terms computed
c     nes     number of solvation energy terms computed
c     nelf    number of metal ligand field energy terms computed
c     neg     number of geometric restraint energy terms computed
c     nex     number of extra energy terms computed
c
c
      module action
      implicit none
      integer neb,nea,neba,neub
      integer neaa,neopb,neopd
      integer neid,neit,net,nept
      integer nebt,neat,nett,nev,nect
      integer nec,necd,ned,nem
      integer nep,new,ner,nes
      integer nelf,neg,nex
      save
      end
