c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module urypot  --  Urey-Bradley functional form details  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     cury       cubic coefficient in Urey-Bradley potential
c     qury       quartic coefficient in Urey-Bradley potential
c     ureyunit   convert Urey-Bradley energy to kcal/mole
c
c
      module urypot
      implicit none
      real*8 cury
      real*8 qury
      real*8 ureyunit
      save
      end
