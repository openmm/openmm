c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module korbs  --  pisystem orbital forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnpi     maximum number of pisystem bond parameter entries
c     maxnpi5    maximum number of 5-membered ring pibond entries
c     maxnpi4    maximum number of 4-membered ring pibond entries
c
c     electron   number of pi-electrons for each atom class
c     ionize     ionization potential for each atom class
c     repulse    repulsion integral value for each atom class
c     sslope     slope for bond stretch vs. pi-bond order
c     sslope5    slope for 5-ring bond stretch vs. pi-bond order
c     sslope4    slope for 4-ring bond stretch vs. pi-bond order
c     tslope     slope for 2-fold torsion vs. pi-bond order
c     tslope5    slope for 5-ring 2-fold torsion vs. pi-bond order
c     tslope4    slope for 4-ring 2-fold torsion vs. pi-bond order
c     kpi        string of atom classes for pisystem bonds
c     kpi5       string of atom classes for 5-ring pisystem bonds
c     kpi4       string of atom classes for 4-ring pisystem bonds
c
c
      module korbs
      use sizes
      implicit none
      integer maxnpi
      integer maxnpi5
      integer maxnpi4
      parameter (maxnpi=500)
      parameter (maxnpi5=200)
      parameter (maxnpi4=200)
      real*8 electron(maxclass)
      real*8 ionize(maxclass)
      real*8 repulse(maxclass)
      real*8 sslope(maxnpi)
      real*8 sslope5(maxnpi5)
      real*8 sslope4(maxnpi4)
      real*8 tslope(maxnpi)
      real*8 tslope5(maxnpi5)
      real*8 tslope4(maxnpi4)
      character*8 kpi(maxnpi)
      character*8 kpi5(maxnpi5)
      character*8 kpi4(maxnpi4)
      save
      end
