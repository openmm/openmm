c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module kbonds  --  bond stretching forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnb   maximum number of bond stretch parameter entries
c     maxnb5  maximum number of 5-membered ring bond stretch entries
c     maxnb4  maximum number of 4-membered ring bond stretch entries
c     maxnb3  maximum number of 3-membered ring bond stretch entries
c     maxnel  maximum number of electronegativity bond corrections
c
c     bcon    force constant parameters for harmonic bond stretch
c     bcon5   force constant parameters for 5-ring bond stretch
c     bcon4   force constant parameters for 4-ring bond stretch
c     bcon3   force constant parameters for 3-ring bond stretch
c     blen    bond length parameters for harmonic bond stretch
c     blen5   bond length parameters for 5-ring bond stretch
c     blen4   bond length parameters for 4-ring bond stretch
c     blen3   bond length parameters for 3-ring bond stretch
c     dlen    electronegativity bond length correction parameters
c     kb      string of atom classes for harmonic bond stretch
c     kb5     string of atom classes for 5-ring bond stretch
c     kb4     string of atom classes for 4-ring bond stretch
c     kb3     string of atom classes for 3-ring bond stretch
c     kel     string of atom classes for electronegativity corrections
c
c
      module kbonds
      implicit none
      integer maxnb
      integer maxnb5
      integer maxnb4
      integer maxnb3
      integer maxnel
      parameter (maxnb=2000)
      parameter (maxnb5=500)
      parameter (maxnb4=500)
      parameter (maxnb3=500)
      parameter (maxnel=500)
      real*8 bcon(maxnb)
      real*8 bcon5(maxnb5)
      real*8 bcon4(maxnb4)
      real*8 bcon3(maxnb3)
      real*8 blen(maxnb)
      real*8 blen5(maxnb5)
      real*8 blen4(maxnb4)
      real*8 blen3(maxnb3)
      real*8 dlen(maxnel)
      character*8 kb(maxnb)
      character*8 kb5(maxnb5)
      character*8 kb4(maxnb4)
      character*8 kb3(maxnb3)
      character*12 kel(maxnel)
      save
      end
