c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ktorsn  --  torsional angle forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnt    maximum number of torsional angle parameter entries
c     maxnt5   maximum number of 5-membered ring torsion entries
c     maxnt4   maximum number of 4-membered ring torsion entries
c
c     t1       torsional parameters for standard 1-fold rotation
c     t2       torsional parameters for standard 2-fold rotation
c     t3       torsional parameters for standard 3-fold rotation
c     t4       torsional parameters for standard 4-fold rotation
c     t5       torsional parameters for standard 5-fold rotation
c     t6       torsional parameters for standard 6-fold rotation
c     t15      torsional parameters for 1-fold rotation in 5-ring
c     t25      torsional parameters for 2-fold rotation in 5-ring
c     t35      torsional parameters for 3-fold rotation in 5-ring
c     t45      torsional parameters for 4-fold rotation in 5-ring
c     t55      torsional parameters for 5-fold rotation in 5-ring
c     t65      torsional parameters for 6-fold rotation in 5-ring
c     t14      torsional parameters for 1-fold rotation in 4-ring
c     t24      torsional parameters for 2-fold rotation in 4-ring
c     t34      torsional parameters for 3-fold rotation in 4-ring
c     t44      torsional parameters for 4-fold rotation in 4-ring
c     t54      torsional parameters for 5-fold rotation in 4-ring
c     t64      torsional parameters for 6-fold rotation in 4-ring
c     kt       string of atom classes for torsional angles
c     kt5      string of atom classes for 5-ring torsions
c     kt4      string of atom classes for 4-ring torsions
c
c
      module ktorsn
      implicit none
      integer maxnt
      integer maxnt5
      integer maxnt4
      parameter (maxnt=2000)
      parameter (maxnt5=500)
      parameter (maxnt4=500)
      real*8 t1(2,maxnt)
      real*8 t2(2,maxnt)
      real*8 t3(2,maxnt)
      real*8 t4(2,maxnt)
      real*8 t5(2,maxnt)
      real*8 t6(2,maxnt)
      real*8 t15(2,maxnt5)
      real*8 t25(2,maxnt5)
      real*8 t35(2,maxnt5)
      real*8 t45(2,maxnt5)
      real*8 t55(2,maxnt5)
      real*8 t65(2,maxnt5)
      real*8 t14(2,maxnt4)
      real*8 t24(2,maxnt4)
      real*8 t34(2,maxnt4)
      real*8 t44(2,maxnt4)
      real*8 t54(2,maxnt4)
      real*8 t64(2,maxnt4)
      character*16 kt(maxnt)
      character*16 kt5(maxnt5)
      character*16 kt4(maxnt4)
      save
      end
