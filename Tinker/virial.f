c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module virial  --  components of internal virial tensor  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     vir         total internal virial Cartesian tensor components
c     use_virial  logical flag governing use of virial computation
c
c
      module virial
      implicit none
      real*8 vir(3,3)
      logical use_virial
      save
      end
