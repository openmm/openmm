c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine field  --  get the potential energy functions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "field" sets the force field potential energy functions from
c     a parameter file and modifications specified in a keyfile
c
c
      subroutine field
      use keys
      use potent
      implicit none
      integer i
      character*240 record
c
c
c     set the default values for the active potentials
c
      use_bond = .true.
      use_angle = .true.
      use_strbnd = .true.
      use_urey = .true.
      use_angang = .true.
      use_opbend = .true.
      use_opdist = .true.
      use_improp = .true.
      use_imptor = .true.
      use_tors = .true.
      use_pitors = .true.
      use_strtor = .true.
      use_angtor = .true.
      use_tortor = .true.
      use_vdw = .true.
      use_ct = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
      use_rxnfld = .false.
      use_solv = .true.
      use_metal = .false.
      use_geom = .true.
      use_extra = .true.
c
c     read the potential energy force field parameter file
c
      call getprm
c
c     check keywords for potential function control parameters
c
      do i = 1, nkey
         record = keyline(i)
         call prmkey (record)
      end do
      return
      end
