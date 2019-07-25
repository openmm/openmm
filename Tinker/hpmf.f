c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module hpmf  --  hydrophobic potential of mean force term  ##
c     ##                                                             ##
c     #################################################################
c
c
c     rcarbon    radius of a carbon atom for use with HPMF
c     rwater     radius of a water molecule for use with HPMF
c     acsurf     surface area of a hydrophobic carbon atom
c     safact     constant for calculation of atomic surface area
c     tslope     tanh slope (set very steep, default=100)
c     toffset    shift the tanh plot along the x-axis (default=6)
c     hpmfcut    cutoff distance for pairwise HPMF interactions
c     h1,h2,h3   hydrophobic PMF well depth parameter
c     c1,c2,c3   hydrophobic PMF well center point
c     w1,w2,w3   reciprocal of the hydrophobic PMF well width
c
c     npmf       number of hydrophobic carbon atoms in the system
c     ipmf       number of the atom for each HPMF carbon atom site
c     rpmf       radius of each atom for use with hydrophobic PMF
c     acsa       SASA value for each hydrophobic PMF carbon atom
c
c
      module hpmf
      implicit none
      real*8 rcarbon,rwater
      real*8 acsurf,safact
      real*8 tslope,toffset
      real*8 hpmfcut
      real*8 h1,h2,h3
      real*8 c1,c2,c3
      real*8 w1,w2,w3
      parameter (rcarbon=1.7d0)
      parameter (rwater=1.4d0)
      parameter (acsurf=120.7628d0)
      parameter (safact=0.3516d0)
      parameter (tslope=100.0d0)
      parameter (toffset=6.0d0)
      parameter (hpmfcut=11.0d0)
      parameter (h1=-0.7308004860404441194d0)
      parameter (h2=0.2001645051578760659d0)
      parameter (h3=-0.0905499953418473502d0)
      parameter (c1=3.8167879266271396155d0)
      parameter (c2=5.4669162286016419472d0)
      parameter (c3=7.1167694861385353278d0)
      parameter (w1=1.6858993102248638341d0)
      parameter (w2=1.3906405621629980285d0)
      parameter (w3=1.5741657341338335385d0)
      integer npmf
      integer, allocatable :: ipmf(:)
      real*8, allocatable :: rpmf(:)
      real*8, allocatable :: acsa(:)
      save
      end
