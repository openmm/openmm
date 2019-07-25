c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module rxnfld  --  reaction field matrix and indices  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ijk   indices into the reaction field element arrays
c     b1    first reaction field matrix element array
c     b2    second reaction field matrix element array
c
c
      module rxnfld
      implicit none
      integer ijk(0:5,0:5,0:5)
      real*8 b1(40,13)
      real*8 b2(40,13)
      save
      end
