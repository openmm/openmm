c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module nonpol  --  nonpolar cavity & dispersion parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     epso      water oxygen eps for implicit dispersion term
c     epsh      water hydrogen eps for implicit dispersion term
c     rmino     water oxygen Rmin for implicit dispersion term
c     rminh     water hydrogen Rmin for implicit dispersion term
c     awater    water number density at standard temp & pressure
c     slevy     enthalpy-to-free energy scale factor for dispersion
c
c     solvprs   limiting microscopic solvent pressure value
c     surften   limiting macroscopic surface tension value
c     spcut     starting radius for solvent pressure tapering
c     spoff     cutoff radius for solvent pressure tapering
c     stcut     starting radius for surface tension tapering
c     stoff     cutoff radius for surface tension tapering
c     rcav      atomic radius of each atom for cavitation energy
c     rdisp     atomic radius of each atom for dispersion energy
c     cdisp     maximum dispersion energy for each atom
c
c
      module nonpol
      implicit none
      real*8 epso,epsh
      real*8 rmino,rminh
      real*8 awater,slevy
      parameter (epso=0.1100d0)
      parameter (epsh=0.0135d0)
      parameter (rmino=1.7025d0)
      parameter (rminh=1.3275d0)
      parameter (awater=0.033428d0)
      parameter (slevy=1.0d0)
      real*8 solvprs,surften
      real*8 spcut,spoff
      real*8 stcut,stoff
      real*8, allocatable :: rcav(:)
      real*8, allocatable :: rdisp(:)
      real*8, allocatable :: cdisp(:)
      save
      end
