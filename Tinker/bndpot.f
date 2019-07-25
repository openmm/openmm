c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module bndpot  --  bond stretch functional form details  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     cbnd      cubic coefficient in bond stretch potential
c     qbnd      quartic coefficient in bond stretch potential
c     bndunit   convert bond stretch energy to kcal/mole
c     bndtyp    type of bond stretch potential energy function
c
c
      module bndpot
      implicit none
      real*8 cbnd
      real*8 qbnd
      real*8 bndunit
      character*8 bndtyp
      save
      end
