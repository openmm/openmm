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

        return energies, tinker, system

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
        energies, _, _ = self.computeAmoebaEnergies(xyzFile, keyFiles)

        # Compare to values computed with Tinker.
        self.assertEnergyEqual(1104.0455, energies["AmoebaBondForce"])
        self.assertEnergyEqual(602.7082, energies["AmoebaAngleForce"] + energies["AmoebaInPlaneAngleForce"])
        self.assertEnergyEqual(2.0572, energies["AmoebaOutOfPlaneBendForce"], 1e-4)
        self.assertEnergyEqual(-0.1361, energies["AmoebaStretchBendForce"], 1e-3)
        self.assertEnergyEqual(-0.8625, energies["PeriodicTorsionForce"], 1e-4)
        self.assertEnergyEqual(-33.8595, energies["HarmonicBondForce"])
        self.assertEnergyEqual(5908.1343, energies["AmoebaVdwForce"])
        self.assertEnergyEqual(-13227.0088 - 5439.9969, energies["AmoebaMultipoleForce"])
        self.assertEnergyEqual(-11084.9187, sum(list(energies.values())))

    def test_Amoeba18Peptide(self):
        """
        Test that TinkerFiles generates a system that gives the same energies as Tinker for a peptide.


        Notes
        -----
        $ analyze systems/peptide.xyz systems/amoebabio18.prm E

        Total Potential Energy :                985.2453 Kcal/mole

        Energy Component Breakdown :           Kcal/mole        Interactions

        Bond Stretching                          19.6519              333
        Angle Bending                            58.2509              596
        Stretch-Bend                             -0.4384              533
        Out-of-Plane Bend                         1.9697              225
        Torsional Angle                          -2.3514              875
        Pi-Orbital Torsion                        1.2115               48
        Torsion-Torsion                          -3.2958                1
        Van der Waals                          1509.1915            52699
        Atomic Multipoles                      -488.0403            52699
        Polarization                           -110.9042            52699

        $ analyze systems/peptide.xyz systems/amoebabio18.prm D
         Overall System Contents :

        Number of Atoms                              328
        Number of Molecules                            1
        Total System Mass                      2396.7620

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
        xyzFile = "systems/peptide.xyz"
        keyFiles = ["systems/amoebabio18.prm"]
        energies, tinker, _ = self.computeAmoebaEnergies(xyzFile, keyFiles)

        # Assert residues are correct
        residues = [residue.name for residue in tinker.topology.residues()]
        assert residues == ['ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO', 'SER', 'THR', 'CYS', 
                            'PHE', 'TYR', 'HIS', 'TRP', 'ASP', 'ASN', 'GLU', 'GLN', 'MET', 
                            'LYS', 'ARG'], f'Unexpected residues: {residues}'

        # Compare to values computed with Tinker.
        self.assertEnergyEqual(19.6519, energies["AmoebaBondForce"])
        self.assertEnergyEqual(58.2509, energies["AmoebaAngleForce"] + energies["AmoebaInPlaneAngleForce"])
        self.assertEnergyEqual(1.9697, energies["AmoebaOutOfPlaneBendForce"], 1e-4)
        self.assertEnergyEqual(-0.4384, energies["AmoebaStretchBendForce"], 1e-3)
        self.assertEnergyEqual(-2.3514, energies["PeriodicTorsionForce"], 1e-4)
        self.assertEnergyEqual(1.2115, energies["AmoebaPiTorsionForce"], 1e-4)
        self.assertEnergyEqual(-3.2958, energies["AmoebaTorsionTorsionForce"])
        self.assertEnergyEqual(1509.1915, energies["AmoebaVdwForce"])
        self.assertEnergyEqual(-488.0403 - 110.9042, energies["AmoebaMultipoleForce"], 1e-3)
        self.assertEnergyEqual(985.2453, sum(list(energies.values())), 1e-3)

    def test_Amoeba18Nucleic(self):
        """
        Test that TinkerFiles generates a system that gives the same energies as Tinker for DNA and RNA.

        Notes
        -----
        $ analyze systems/nucleic.xyz systems/amoebabio18.prm E

        Intermolecular Energy :                 896.3435 Kcal/mole

        Total Potential Energy :               3146.3045 Kcal/mole

        Energy Component Breakdown :           Kcal/mole        Interactions

        Bond Stretching                         749.6953              827
        Angle Bending                           579.9971             1483
        Stretch-Bend                              5.2225             1253
        Out-of-Plane Bend                        10.6630              441
        Torsional Angle                         166.7233             2197
        Pi-Orbital Torsion                       57.2066              142
        Stretch-Torsion                          -4.2538               68
        Angle-Torsion                            -5.0402              112
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
        energies, tinker, _ = self.computeAmoebaEnergies(xyzFile, keyFiles)

        # Assert residues are correct
        residues = [residue.name for residue in tinker.topology.residues()]
        assert residues == ['DA', 'DG', 'DC', 'DG', 'DT', 'DG', 'DG', 'DG', 'DA', 'DC', 'DC', 
                            'G', 'C', 'G', 'U', 'U', 'A', 'A', 'G', 'U', 'C', 'G', 'C', 'A'], f'Unexpected residues: {residues}'

        # Compare to values computed with Tinker.
        self.assertEnergyEqual(749.6953, energies["AmoebaBondForce"])
        self.assertEnergyEqual(579.9971, energies["AmoebaAngleForce"] + energies["AmoebaInPlaneAngleForce"])
        self.assertEnergyEqual(10.6630, energies["AmoebaOutOfPlaneBendForce"])
        self.assertEnergyEqual(5.2225, energies["AmoebaStretchBendForce"])
        self.assertEnergyEqual(166.7233, energies["PeriodicTorsionForce"])
        self.assertEnergyEqual(57.2066, energies["AmoebaPiTorsionForce"])
        self.assertEnergyEqual(-4.2538, energies["AmoebaStretchTorsionForce"])
        self.assertEnergyEqual(-5.0402, energies["AmoebaAngleTorsionForce"], 1e-3)
        self.assertEnergyEqual(187.1103, energies["AmoebaVdwForce"])
        self.assertEnergyEqual(1635.1289 - 236.1484, energies["AmoebaMultipoleForce"])
        self.assertEnergyEqual(3146.3045, sum(list(energies.values())))

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

    def test_Topology(self):
        """Test that the Topology created from Tinker files is correct."""
        xyzFile = "systems/ubiquitin.xyz"
        keyFiles = ["systems/amoebabio18.prm"]
        tinker = TinkerFiles(xyzFile, keyFiles)
        topology = tinker.topology

        # Test ubiquitin 
        assert topology.getNumAtoms() == 1406, f'Expected 1406 atoms for ubiquitin, found {topology.getNumAtoms()}'
        residues = [residue.name for residue in topology.residues()]
        assert residues == ["MET", "GLN", "ILE", "PHE", "VAL", "LYS", "THR", "LEU", "THR", "GLY", "LYS", "THR", "ILE", "THR", "LEU",
                            "GLU", "VAL", "GLU", "PRO", "SER", "ASP", "THR", "ILE", "GLU", "ASN", "VAL", "LYS", "ALA", "LYS", "ILE",
                            "GLN", "ASP", "LYS", "GLU", "GLY", "ILE", "PRO", "PRO", "ASP", "GLN", "GLN", "ARG", "LEU", "ILE", "PHE",
                            "ALA", "GLY", "LYS", "GLN", "LEU", "GLU", "ASP", "GLY", "ARG", "THR", "LEU", "SER", "ASP", "TYR", "ASN",
                            "ILE", "GLN", "LYS", "GLU", "SER", "THR", "LEU", "HIS", "LEU", "VAL", "LEU", "ARG", "LEU", "ARG", "GLY", "GLY"] + ["HOH"]*58, f'Unexpected residues: {residues}'
            
        xyzFile = "systems/ubiquitin.xyz"
        keyFiles = ["systems/amoebabio18.prm"]
        tinker = TinkerFiles(xyzFile, keyFiles)
        topology = tinker.topology

        # Test bdna
        xyzFile = "systems/bdna.xyz"

        tinker = TinkerFiles(xyzFile, keyFiles)
        topology = tinker.topology
        assert topology.getNumAtoms() == 758, f'Expected 758 atoms for bdna, found {topology.getNumAtoms()}'
        residues = [residue.name for residue in topology.residues()]
        assert residues == ["DC", "DG", "DC", "DG", "DA", "DA", "DT", "DT", "DC", "DG", "DC", "DG"] * 2, f'Unexpected residues: {residues}' 
        assert topology.getNumChains() == 2, f'Expected 2 chains for bdna, found {topology.getNumChains()}'

   