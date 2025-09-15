import unittest

from openmm import *
from openmm.app import *
from openmm.unit import *


class TestTinkerFiles(unittest.TestCase):
    """Test the TinkerFiles class for reading Tinker files and creating OpenMM systems."""

    def computeAmoebaEnergies(self, xyzFile, keyFiles):
        tinker = TinkerFiles(xyzFile, keyFiles)
        system = tinker.createSystem(
            polarization="mutual",
            mutualInducedTargetEpsilon=1e-5,
            nonbondedMethod=NoCutoff,
            constraints=None,
            useDispersionCorrection=False,
        )

        # Compute the energy with OpenMM.
        for i, f in enumerate(system.getForces()):
            f.setForceGroup(i)

        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform("Reference"))
        context.setPositions(tinker.getPositions())
        energies = {}
        for i, f in enumerate(system.getForces()):
            state = context.getState(getEnergy=True, groups={i})
            energies[f.getName()] = state.getPotentialEnergy().value_in_unit(
                kilocalories_per_mole
            )

        return energies

    def assertEnergyEqual(self, expected, found, rtol=1e-5):
        """Assert the values match to a specified relative precision."""
        if abs(expected-found) > rtol*abs(expected):
            raise AssertionError(f'{expected} != {found} to a relative tolerance of {rtol}')

    def test_Amoeba18PoltypePhenolWater(self):
        """
        Test that TinkerFiles generates a system that gives the same energies as Tinker for phenol (poltype) in water (amoebabio18).

        Notes
        -----
        $ analyze systems/phenol_water.xyz systems/test.prm E
        (test.prm is $ cat systems/phenol.prm systems/amoebabio18.prm > systems/test.prm)

        Intermolecular Energy :              -12759.7770 Kcal/mole

        Total Potential Energy :             -11084.9187 Kcal/mole

        Energy Component Breakdown :           Kcal/mole        Interactions

        Bond Stretching                        1104.0455             3007
        Angle Bending                           602.7082             1516
        Stretch-Bend                             -0.1361               19
        Urey-Bradley                            -33.8595             1497
        Out-of-Plane Bend                         2.0572               18
        Torsional Angle                          -0.8625               26
        Van der Waals                          5908.1343         10136233
        Atomic Multipoles                    -13227.0088         10136233
        Polarization                          -5439.9969         10136233

        $ analyze systems/phenol_water.xyz systems/test.prm D
        (test.prm is $ cat systems/phenol.prm systems/amoebabio18.prm > systems/test.prm)

        Overall System Contents :

        Number of Atoms                             4504
        Number of Molecules                         1498
        Total System Mass                     27062.5680

        Force Field Name :               AMOEBA-BIO-2018

        VDW Function                       BUFFERED-14-7
        Size Descriptor                            R-MIN
        Size Unit Type                          DIAMETER
        Size Combining Rule                   CUBIC-MEAN
        Well Depth Rule                              HHG

        Electrostatics                  ATOMIC MULTIPOLE
        Induction                         INDUCED DIPOLE
        Polarization                              MUTUAL
        Polarization                       THOLE DAMPING
        """
        xyzFile = "systems/phenol_water.xyz"
        keyFiles = ["systems/phenol.prm", "systems/amoebabio18.prm"]
        energies = self.computeAmoebaEnergies(xyzFile, keyFiles)

        # Compare to values computed with Tinker.
        self.assertEnergyEqual(1104.0455, energies["AmoebaBond"])
        self.assertEnergyEqual(602.7082, energies["AmoebaAngle"] + energies["AmoebaInPlaneAngle"])
        self.assertEnergyEqual(2.0572, energies["AmoebaOutOfPlaneBend"])
        self.assertEnergyEqual(-0.1361, energies["AmoebaStretchBend"])
        self.assertEnergyEqual(-0.8625, energies["PeriodicTorsionForce"])
        self.assertEnergyEqual(-33.8595, energies["HarmonicBondForce"])
        self.assertEnergyEqual(5908.1343, energies["AmoebaVdwForce"])
        self.assertEnergyEqual(-13227.0088 - 5439.9969, energies["AmoebaMultipoleForce"])
        self.assertEnergyEqual(-11084.9187, sum(list(energies.values())))

    def test_Amoeba18Nucleic(self):
        """
        Test that TinkerFiles generates a system that gives the same energies as Tinker for DNA and RNA.

        Notes
        -----
        $ analyze systems/nucleic.xyz systems/amoebabio18.prm D

        Intermolecular Energy :                 896.3435 Kcal/mole

        Total Potential Energy :               3151.2568 Kcal/mole

        Energy Component Breakdown :           Kcal/mole        Interactions

        Bond Stretching                         749.6953              827
        Angle Bending                           579.9971             1483
        Stretch-Bend                              5.2225             1253
        Out-of-Plane Bend                        10.6630              441
        Torsional Angle                         166.7233             2197
        Pi-Orbital Torsion                       57.2066              142
        Stretch-Torsion                          -4.2538               68
        Angle-Torsion                            -0.0880              112
        Van der Waals                           187.1103           291451
        Atomic Multipoles                      1635.1289           291451
        Polarization                           -236.1484           291451

        $ analyze systems/nucleic.xyz systems/amoebabio18.prm G

        Overall System Contents :

        Number of Atoms                              767
        Number of Molecules                            2
        Total System Mass                      7500.6170

        Force Field Name :               AMOEBA-BIO-2018

        VDW Function                       BUFFERED-14-7
        Size Descriptor                            R-MIN
        Size Unit Type                          DIAMETER
        Size Combining Rule                   CUBIC-MEAN
        Well Depth Rule                              HHG

        Electrostatics                  ATOMIC MULTIPOLE
        Induction                         INDUCED DIPOLE
        Polarization                              MUTUAL
        Polarization                       THOLE DAMPING
        """
        xyzFile = "systems/nucleic.xyz"
        keyFiles = ["systems/amoebabio18.prm"]
        energies = self.computeAmoebaEnergies(xyzFile, keyFiles)

        # Compare to values computed with Tinker.
        self.assertEnergyEqual(749.6953, energies["AmoebaBond"])
        self.assertEnergyEqual(579.9971, energies["AmoebaAngle"] + energies["AmoebaInPlaneAngle"])
        self.assertEnergyEqual(10.6630, energies["AmoebaOutOfPlaneBend"])
        self.assertEnergyEqual(5.2225, energies["AmoebaStretchBend"])
        self.assertEnergyEqual(166.7233, energies["PeriodicTorsionForce"])
        self.assertEnergyEqual(57.2066, energies["AmoebaPiTorsion"])
        self.assertEnergyEqual(-4.2538, energies["AmoebaStretchTorsion"])
        self.assertEnergyEqual(-0.0880, energies["AmoebaAngleTorsion"], 1e-4)
        self.assertEnergyEqual(187.1103, energies["AmoebaVdwForce"])
        self.assertEnergyEqual(1635.1289 - 236.1484, energies["AmoebaMultipoleForce"])
        self.assertEnergyEqual(3151.2568, sum(list(energies.values())))

    def test_Amoeba13ForcesImplicit(self):
        """Compute forces and compare them to ones generated with a previous version of OpenMM to ensure they haven't changed."""
        # Define mapping between positions of the .xyz file created using
        # pdbxyz alanine-dipeptide-implicit.pdb
        # and the original .pdb file.
        mapping = {1: 4, 2: 1, 3: 5, 4: 6, 5: 2, 6: 3, 7: 7, 8: 11, 
                   9: 8, 10: 12, 11: 13, 12: 14, 13: 15, 14: 16, 15: 9, 
                   16: 10, 17: 17, 18: 19, 19: 18, 20: 20, 21: 21, 22: 22 
        }

        xyzFile = 'systems/alanine-dipeptide-implicit.xyz'
        keyFiles = ['systems/amoebapro13.prm']

        tinker = TinkerFiles(xyzFile, keyFiles)
        system = tinker.createSystem(polarization='direct',
                                     constraints=None, 
                                     implicitSolvent=True,
        )

        # Assert residues are correct
        residues = [residue.name for residue in tinker.topology.residues()]
        assert residues == ['ACE', 'ALA', 'NME'], f'Unexpected residues: {residues}'

        # Compute the forces with OpenMM
        integrator = VerletIntegrator(0.001)
        context = Context(system, integrator, Platform.getPlatform('Reference'))
        context.setPositions(tinker.getPositions())
        state1 = context.getState(getForces=True)

        with open('systems/alanine-dipeptide-amoeba-forces.xml') as input:
            state2 = XmlSerializer.deserialize(input.read())

        state1Forces = state1.getForces().value_in_unit(kilojoules_per_mole/nanometer)
        state2Forces = state2.getForces().value_in_unit(kilojoules_per_mole/nanometer)

        for i, j in mapping.items():
            f1 = state1Forces[j-1]
            f2 = state2Forces[i-1]
            diff = norm(f1-f2)
            self.assertTrue(diff < 0.1 or diff/norm(f1) < 1e-3, f1-f2)
