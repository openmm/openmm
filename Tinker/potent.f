c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module potent  --  usage of potential energy components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     use_bond    logical flag governing use of bond stretch potential
c     use_angle   logical flag governing use of angle bend potential
c     use_strbnd  logical flag governing use of stretch-bend potential
c     use_urey    logical flag governing use of Urey-Bradley potential
c     use_angang  logical flag governing use of angle-angle cross term
c     use_opbend  logical flag governing use of out-of-plane bend term
c     use_opdist  logical flag governing use of out-of-plane distance
c     use_improp  logical flag governing use of improper dihedral term
c     use_imptor  logical flag governing use of improper torsion term
c     use_tors    logical flag governing use of torsional potential
c     use_pitors  logical flag governing use of pi-system torsion term
c     use_strtor  logical flag governing use of stretch-torsion term
c     use_angtor  logical flag governing use of angle-torsion term
c     use_tortor  logical flag governing use of torsion-torsion term
c     use_vdw     logical flag governing use of vdw der Waals potential
c     use_ct      logical flag governing use of charge transfer potential
c     use_charge  logical flag governing use of charge-charge potential
c     use_chgdpl  logical flag governing use of charge-dipole potential
c     use_dipole  logical flag governing use of dipole-dipole potential
c     use_mpole   logical flag governing use of multipole potential
c     use_polar   logical flag governing use of polarization term
c     use_rxnfld  logical flag governing use of reaction field term
c     use_solv    logical flag governing use of continuum solvation
c     use_metal   logical flag governing use of ligand field term
c     use_geom    logical flag governing use of geometric restraints
c     use_extra   logical flag governing use of extra potential term
c     use_born    logical flag governing use of Born radii values
c     use_orbit   logical flag governing use of pisystem computation
c
c
      module potent
      implicit none
      logical use_bond,use_angle
      logical use_strbnd,use_urey
      logical use_angang,use_opbend
      logical use_opdist,use_improp
      logical use_imptor,use_tors
      logical use_pitors,use_strtor
      logical use_angtor,use_tortor
      logical use_vdw,use_charge
      logical use_ct
      logical use_chgdpl,use_dipole
      logical use_mpole,use_polar
      logical use_rxnfld,use_solv
      logical use_metal,use_geom
      logical use_extra,use_born
      logical use_orbit
      save
      end
