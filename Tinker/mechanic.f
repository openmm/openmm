c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic
      use inform
      use iounit
      use limits
      use potent
      use vdwpot
      use ctran
      implicit none
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign bond, angle and cross term potential parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor .or.
     &    (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or.use_angang .or. 
     &    use_angtor .or. use_opbend .or. use_opdist)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend
      if (use_angle .or. use_opdist)  call kopdist
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or.
     &    use_angtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors
      if (use_strtor)  call kstrtor
      if (use_angtor)  call kangtor
      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
      if (use_ct)  call kctr
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld)  call kmpole
      if (use_mpole .or. use_polar .or. use_solv)  call kpolar
      if (use_ewald)  call kewald
c
c     assign solvation, metal, pisystem and restraint parameters
c
      if (use_solv)  call ksolv
      if (use_metal)  call kmetal
      if (use_orbit)  call korbit
      if (use_geom)  call kgeom
      if (use_extra)  call kextra
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
