.. _amoeba-plugin:

AMOEBA Plugin
#############

OpenMM |version| provides a plugin that implements the AMOEBA polarizable atomic
multipole force field from Jay Ponder’s lab. The AMOEBA force field may be used
through OpenMM’s Python application layer. We have also created a modified
version of TINKER (referred to as TINKER-OpenMM here) that uses OpenMM to
accelerate AMOEBA simulations. TINKER-OpenMM can be created from a TINKER
package using three files made available through the OpenMM home page. OpenMM
AMOEBA Force and System objects containing AMOEBA forces can be serialized.

In the following sections, the individual forces and options available in the
plugin are listed, and the steps required to build and use the plugin and
TINKER-OpenMM are outlined. Validation results are also reported.  Benchmarks
can be found on the OpenMM wiki at http://wiki.simtk.org/openmm/Benchmarks.

OpenMM AMOEBA Supported Forces and Options
*******************************************


.. _supported-forces-and-options:

Supported Forces and Options
============================

The AMOEBA force terms implemented in OpenMM are listed in :autonumref:`Table,mapping from TINKER` along
with the supported and unsupported options. TINKER options that are not
supported for any OpenMM force include the grouping of atoms (e.g. protein
chains), the infinite polymer check, and no exclusion of particles from
energy/force calculations (‘active’/’inactive’ particles).  The virial is not
calculated for any force.

All rotation axis types are supported: ‘Z-then-X’, ‘Bisector’, ‘Z-Bisect’,
‘3-Fold’, ‘Z-Only’.


=================================  ==================================  ======================================================================================================================================================================================
TINKER Force                       OpenMM Force                        Option/Note
=================================  ==================================  ======================================================================================================================================================================================
ebond1 (bondterm)                  AmoebaBondForce                     bndtyp='HARMONIC' supported, 'MORSE' not implemented
Eangle71 (angleterm)               AmoebaAngleForce                    angtyp='HARMONIC' and 'IN-PLANE' supported; 'LINEAR' and 'FOURIER' not implemented
etors1a (torsionterm)              PeriodicTorsionForce                All options implemented; smoothing version(etors1b) not supported
etortor1 (tortorterm)              AmoebaTorsionTorsionForce           All options implemented
eopbend1 (opbendterm)              AmoebaOutOfPlaneBendForce           opbtyp = 'ALLINGER' implemented; 'W-D-C' not implemented
epitors1 (pitorsterm)              AmoebaPiTorsionForce                All options implemented
estrbnd1 (strbndterm)              AmoebaStretchBendForce              All options implemented
ehal1a (vdwterm)                   AmoebaVdwForce                      ehal1b(LIGHTS) not supported
empole1a (mpoleterm)               AmoebaMultipoleForce                poltyp = 'MUTUAL', 'DIRECT'  supported
empole1c (mpoleterm) PME           AmoebaMultipoleForce                poltyp = 'MUTUAL', 'DIRECT' supported; boundary= 'VACUUM' unsupported
esolv1 (solvateterm)               | AmoebaWcaDispersionForce,         Only born-radius=’grycuk’ and solvate=’GK’ supported; unsupported solvate settings:
                                   | AmoebaGeneralizedKirkwoodForce    ‘ASP’, ‘SASA’, ‘ONION’, ‘pb’, 'GB-HPMF’, 'Gk-HPMF’; SASA computation is based on ACE approximation
eurey1 (ureyterm)                  HarmonicBondForce                   All options implemented
=================================  ==================================  ======================================================================================================================================================================================

:autonumber:`Table,mapping from TINKER`\ :  Mapping between TINKER and OpenMM AMOEBA forces


Some specific details to be aware of are the following:

* Forces available in TINKER but not implemented in the OpenMM AMOEBA plugin
  include the following: angle-angle, out-of-plane distance, improper dihedral,
  improper torsion, stretch-torsion, charge-charge, atomwise charge-dipole,
  dipole-dipole, reaction field, ligand field, restraint, scf molecular orbital
  calculation; strictly speaking, these are not part of the AMOEBA force field.

* Implicit solvent in TINKER-OpenMM is implemented with key file entry ‘solvate
  GK’.  The entry ‘born-radius grycuk’ should also be included; only the ‘grycuk’
  option for calculating the Born radii is available in the plugin.

* In TINKER, the nonpolar cavity contribution to the solvation term is
  calculated using an algorithm that does not map well to GPUs.  Instead the
  OpenMM plugin uses the TINKER version of the ACE approximation to estimate the
  cavity contribution to the SASA.

* Calculations using the CUDA platform may be done in either single or double
  precision; for the Reference platform, double precision is used.  TINKER uses
  double precision.

* The TINKER parameter files for the AMOEBA force-field parameters are based on
  units of kilocalorie/Å, whereas OpenMM uses units of kilojoules/nanometer; both
  TINKER and OpenMM use picoseconds time units. Hence, in mapping the force-field
  parameters from TINKER files to OpenMM, many of the parameter values must be
  converted to the OpenMM units. The setup methods in the TINKER-OpenMM
  application perform the required conversions.


Supported Integrators
=====================

In addition to the limitations to the forces outlined above, TINKER-OpenMM can
only use either the ‘Verlet’ or ‘Stochastic’ integrators when the OpenMM plugin
is used; an equivalent to the TINKER ‘Beeman’ integrator is unavailable in
OpenMM.

OpenMM AMOEBA Validation
************************

OpenMM and TINKER 6.1.01 were each used to compute the atomic forces for
dihydrofolate reductase (DHFR) in implicit and explicit solvent.  Calculations
used the CUDA platform, and were repeated for both single and double precision.
For every atom, the relative difference between OpenMM and TINKER was computed
as 2·\|F\ :sub:`MM`\ –F\ :sub:`T`\ \|/(\|F\ :sub:`MM`\ \|+\|F\ :sub:`T`\ \|), where
F\ :sub:`MM` is the force computed by OpenMM and F\ :sub:`T` is the force
computed by TINKER.  The median over all atoms is shown in :autonumref:`Table,comparison to TINKER`\ .

Because OpenMM and TINKER use different approximations to compute the cavity
term, the differences in forces are much larger for implicit solvent than for
explicit solvent.  We therefore repeated the calculations, removing the cavity
term.  This yields much closer agreement between OpenMM and TINKER,
demonstrating that the difference comes entirely from that one term.

=========================  ==========================  ===================
Solvent Model              single                      double
=========================  ==========================  ===================
Implicit                   1.04·10\ :sup:`-2`          1.04·10\ :sup:`-2`
Implicit (no cavity term)  9.23·10\ :sup:`-6`          1.17·10\ :sup:`-6`
Explicit                   3.73·10\ :sup:`-5`          1.83·10\ :sup:`-7`
=========================  ==========================  ===================

:autonumber:`Table,comparison to TINKER`\ :  Median relative difference in forces between OpenMM and TINKER

