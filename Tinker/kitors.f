c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kitors  --  improper torsion forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnti   maximum number of improper torsion parameter entries
c
c     ti1      torsional parameters for improper 1-fold rotation
c     ti2      torsional parameters for improper 2-fold rotation
c     ti3      torsional parameters for improper 3-fold rotation
c     kti      string of atom classes for improper torsional parameters
c
c
      module kitors
      implicit none
      integer maxnti
      parameter (maxnti=500)
      real*8 ti1(2,maxnti)
      real*8 ti2(2,maxnti)
      real*8 ti3(2,maxnti)
      character*16 kti(maxnti)
      save
      end
