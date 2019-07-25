c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ksttor  --  stretch-torsion forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnbt   maximum number of stretch-torsion parameter entries
c
c     btcon    torsional amplitude parameters for stretch-torsion
c     kbt      string of atom classes for stretch-torsion terms
c
c
      module ksttor
      implicit none
      integer maxnbt
      parameter (maxnbt=500)
      real*8 btcon(9,maxnbt)
      character*16 kbt(maxnbt)
      save
      end
