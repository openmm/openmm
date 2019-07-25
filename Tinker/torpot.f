c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module torpot  --  torsional functional form details  ##
c     ##                                                        ##
c     ############################################################
c
c
c     idihunit  convert improper dihedral energy to kcal/mole
c     itorunit  convert improper torsion amplitudes to kcal/mole
c     torsunit  convert torsional parameter amplitudes to kcal/mole
c     ptorunit  convert pi-system torsion energy to kcal/mole
c     storunit  convert stretch-torsion energy to kcal/mole
c     atorunit  convert angle-torsion energy to kcal/mole
c     ttorunit  convert torsion-torsion energy to kcal/mole
c
c
      module torpot
      implicit none
      real*8 idihunit
      real*8 itorunit
      real*8 torsunit
      real*8 ptorunit
      real*8 storunit
      real*8 atorunit
      real*8 ttorunit
      save
      end
