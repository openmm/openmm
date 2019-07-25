c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2014 by Alex Albaugh & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module ielscf  --  extended Lagrangian induced dipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nfree_aux    total degrees of freedom for auxiliary dipoles
c     tautemp_aux  time constant for auliliary Berendsen thermostat
c     kelvin_aux   target system temperature for auxiliary dipoles
c     uaux         auxiliary induced dipole value at each site
c     upaux        auxiliary shadow induced dipoles at each site
c     vaux         auxiliary induced dipole velocity at each site
c     vpaux        auxiliary shadow dipole velocity at each site
c     aaux         auxiliary induced dipole acceleration at each site
c     apaux        auxiliary shadow dipole acceleration at each site
c     use_ielscf   flag to use inertial extended Lagrangian method
c
c
      module ielscf
      implicit none
      integer nfree_aux
      real*8 tautemp_aux
      real*8 kelvin_aux
      real*8, allocatable :: uaux(:,:)
      real*8, allocatable :: upaux(:,:)
      real*8, allocatable :: vaux(:,:)
      real*8, allocatable :: vpaux(:,:)
      real*8, allocatable :: aaux(:,:)
      real*8, allocatable :: apaux(:,:)
      logical use_ielscf
      save
      end
