c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module energi  --  individual potential energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     esum   total potential energy of the system
c     eb     bond stretch potential energy of the system
c     ea     angle bend potential energy of the system
c     eba    stretch-bend potential energy of the system
c     eub    Urey-Bradley potential energy of the system
c     eaa    angle-angle potential energy of the system
c     eopb   out-of-plane bend potential energy of the system
c     eopd   out-of-plane distance potential energy of the system
c     eid    improper dihedral potential energy of the system
c     eit    improper torsion potential energy of the system
c     et     torsional potential energy of the system
c     ept    pi-system torsion potential energy of the system
c     ebt    stretch-torsion potential energy of the system
c     eat    angle-torsion potential energy of the system
c     ett    torsion-torsion potential energy of the system
c     ev     van der Waals potential energy of the system
CW!!!
c     ect    charge transfer potential energy of the system
CW!!!END

c     ec     charge-charge potential energy of the system
c     ecd    charge-dipole potential energy of the system
c     ed     dipole-dipole potential energy of the system
c     em     atomic multipole potential energy of the system
c     ep     polarization potential energy of the system
c     er     reaction field potential energy of the system
c     es     solvation potential energy of the system
c     elf    metal ligand field potential energy of the system
c     eg     geometric restraint potential energy of the system
c     ex     extra term potential energy of the system
c
c
      module energi
      implicit none
      real*8 esum,eb,ea
      real*8 eba,eub,eaa
      real*8 eopb,eopd,eid
      real*8 eit,et,ept
      real*8 ebt,eat,ett
      real*8 ev,ec,ecd
CW!!!
      real*8 ect
CW!!!END
      real*8 ed,em,ep
      real*8 er,es,elf
      real*8 eg,ex
      save
      end
