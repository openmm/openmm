c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module solute  --  continuum solvation model parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     doffset   dielectric offset to continuum solvation atomic radii
c     p1        single-atom scale factor for analytical Still radii
c     p2        1-2 interaction scale factor for analytical Still radii
c     p3        1-3 interaction scale factor for analytical Still radii
c     p4        nonbonded scale factor for analytical Still radii
c     p5        soft cutoff parameter for analytical Still radii
c     rsolv     atomic radius of each atom for continuum solvation
c     asolv     atomic surface area solvation parameters
c     rborn     Born radius of each atom for GB/SA solvation
c     drb       solvation derivatives with respect to Born radii
c     drbp      GK polarization derivatives with respect to Born radii
c     drobc     chain rule term for Onufriev-Bashford-Case radii
c     gpol      polarization self-energy values for each atom
c     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii
c     aobc      alpha values for Onufriev-Bashford-Case radii
c     bobc      beta values for Onufriev-Bashford-Case radii
c     gobc      gamma values for Onufriev-Bashford-Case radii
c     vsolv     atomic volume of each atom for use with ACE
c     wace      "omega" values for atom class pairs for use with ACE
c     s2ace     "sigma^2" values for atom class pairs for use with ACE
c     uace      "mu" values for atom class pairs for use with ACE
c     solvtyp   type of continuum solvation energy model in use
c     borntyp   method to be used for the Born radius computation
c
c
      module solute
      implicit none
      real*8 doffset
      real*8 p1,p2,p3,p4,p5
      real*8, allocatable :: rsolv(:)
      real*8, allocatable :: asolv(:)
      real*8, allocatable :: rborn(:)
      real*8, allocatable :: drb(:)
      real*8, allocatable :: drbp(:)
      real*8, allocatable :: drobc(:)
      real*8, allocatable :: gpol(:)
      real*8, allocatable :: shct(:)
      real*8, allocatable :: aobc(:)
      real*8, allocatable :: bobc(:)
      real*8, allocatable :: gobc(:)
      real*8, allocatable :: vsolv(:)
      real*8, allocatable :: wace(:,:)
      real*8, allocatable :: s2ace(:,:)
      real*8, allocatable :: uace(:,:)
      character*8 solvtyp
      character*8 borntyp
      save
      end
