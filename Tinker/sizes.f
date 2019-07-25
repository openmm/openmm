c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module sizes  --  parameters to set array dimensions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "sizes" sets values for critical array dimensions used
c     throughout the software; these parameters fix the size of
c     the largest systems that can be handled
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxval          atoms directly bonded to an atom
c     maxref          stored reference molecular systems
c     maxgrp          user-defined groups of atoms
c     maxres          residues in the macromolecule
c     maxfix          geometric constraints and restraints
c
c
      module sizes
      implicit none
      integer maxatm,maxtyp
      integer maxclass,maxval
      integer maxref,maxgrp
      integer maxres,maxfix
      parameter (maxatm=1000000)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxval=8)
      parameter (maxref=30)
      parameter (maxgrp=1000)
      parameter (maxres=10000)
      parameter (maxfix=100000)
      save
      end
