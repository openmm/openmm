c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module paths  --  Elber reaction path method parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     pnorm    length of the reactant-product vector
c     acoeff   transformation matrix 'A' from Elber algorithm
c     pc0      reactant Cartesian coordinates as variables
c     pc1      product Cartesian coordinates as variables
c     pvect    vector connecting the reactant and product
c     pstep    step per cycle along reactant-product vector
c     pzet     current projection on reactant-product vector
c     gc       gradient of the path constraints
c
c
      module paths
      implicit none
      real*8 pnorm
      real*8 acoeff(7,7)
      real*8, allocatable :: pc0(:)
      real*8, allocatable :: pc1(:)
      real*8, allocatable :: pvect(:)
      real*8, allocatable :: pstep(:)
      real*8, allocatable :: pzet(:)
      real*8, allocatable :: gc(:,:)
      save
      end
