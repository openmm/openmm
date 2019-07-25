c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module fields  --  molecular mechanics force field type  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxbio       maximum number of biopolymer atom definitions
c
c     biotyp       force field atom type of each biopolymer type
c     forcefield   string used to describe the current forcefield
c
c
      module fields
      implicit none
      integer maxbio
      parameter (maxbio=10000)
      integer biotyp(maxbio)
      character*20 forcefield
      save
      end
