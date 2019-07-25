c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module xtals  --  structures used for parameter fitting  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxlsq       maximum number of least squares variables
c     maxrsd       maximum number of residual functions
c
c     nxtal        number of molecular structures to be stored
c     nvary        number of potential parameters to optimize
c     ivary        index for the types of potential parameters
c     iresid       structure to which each residual function refers
c     vary         atom numbers involved in potential parameters
c     e0_lattice   ideal lattice energy for the current crystal
c     vartyp       type of each potential parameter to be optimized
c     rsdtyp       experimental variable for each of the residuals
c
c
      module xtals
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=1000)
      parameter (maxrsd=1000)
      integer nxtal,nvary
      integer ivary(maxlsq)
      integer iresid(maxrsd)
      integer vary(2,maxlsq)
      real*8 e0_lattice
      character*16 vartyp(maxlsq)
      character*16 rsdtyp(maxrsd)
      save
      end
