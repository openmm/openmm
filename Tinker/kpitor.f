c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kpitor  --  pi-system torsion forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnpt   maximum number of pi-system torsion parameter entries
c
c     ptcon    force constant parameters for pi-system torsions
c     kpt      string of atom classes for pi-system torsion terms
c
c
      module kpitor
      implicit none
      integer maxnpt
      parameter (maxnpt=500)
      real*8 ptcon(maxnpt)
      character*8 kpt(maxnpt)
      save
      end
