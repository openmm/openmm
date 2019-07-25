c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module polar  --  induced dipole moments & polarizability  ##        
c     ##                                                             ##
c     #################################################################
c
c
c     maxopt    maximum order for OPT induced dipole extrapolation
c
c     npolar    total number of polarizable sites in the system
c     coptmax   maximum coefficient order for OPT dipole extrapolation
c     optlevel  current OPT order for reciprocal potential and field
c     copt      coefficients for OPT total induced dipole moments
c     copm      coefficients for OPT incremental induced dipole moments
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     thole     Thole polarizability damping value for each site
c     dirdamp   direct polarizability damping value for each site
c     penalpha  penetration damping value for each site
c     pdamp     value of polarizability scale factor for each site
c     udir      direct induced dipole components at each multipole site
c     udirp     direct induced dipoles in field used for energy terms
c     udirs     direct GK or PB induced dipoles at each multipole site
c     udirps    direct induced dipoles in field used for GK or PB energy
c     uind      mutual induced dipole components at each multipole site
c     uinp      mutual induced dipoles in field used for energy terms
c     uinds     mutual GK or PB induced dipoles at each multipole site
c     uinps     mutual induced dipoles in field used for GK or PB energy
c     uopt      OPT induced dipole components at each multipole site
c     uoptp     OPT induced dipoles in field used for energy terms
c     uopts     OPT GK or PB induced dipoles at each multipole site
c     uoptps    OPT induced dipoles in field used for GK or PB energy
c     fopt      OPT fractional reciprocal potentials at multipole sites
c     foptp     OPT fractional reciprocal potentials for energy terms
c     uexact    exact SCF induced dipoles to full numerical precision
c     douind    flag to allow induced dipoles at each atomic site
c

      module polar
      implicit none
      integer maxopt
      parameter (maxopt=4)
      integer npolar
      integer coptmax
      integer optlevel
      real*8, allocatable :: copt(:)
      real*8, allocatable :: dirdamp(:) 
      real*8, allocatable :: penalpha(:) 
      real*8, allocatable :: copm(:)
      real*8, allocatable :: polarity(:)
      real*8, allocatable :: thole(:)
      real*8, allocatable :: pdamp(:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: udirs(:,:)
      real*8, allocatable :: udirps(:,:)
      real*8, allocatable :: uind(:,:)
      real*8, allocatable :: uinp(:,:)
      real*8, allocatable :: uinds(:,:)
      real*8, allocatable :: uinps(:,:)
      real*8, allocatable :: uopt(:,:,:)
      real*8, allocatable :: uoptp(:,:,:)
      real*8, allocatable :: uopts(:,:,:)
      real*8, allocatable :: uoptps(:,:,:)
      real*8, allocatable :: fopt(:,:,:)
      real*8, allocatable :: foptp(:,:,:)
      real*8, allocatable :: uexact(:,:)
      logical, allocatable :: douind(:)
      save
      end
