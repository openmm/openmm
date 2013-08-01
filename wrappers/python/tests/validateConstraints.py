from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def validateConstraints(self, topology, system, constraints_value, rigidWater_value):
    """ Given a Topology, System, a value for 'constraints' and a value for
    'rigidWater', verify that the System contains the correct bonds and constraints.    

    """    
    
    # Build the set of expected constraints.
    expected_constraint_set=set()
    if rigidWater_value:
        for residue in topology.residues():
            if residue.name=="HOH":
                indices = [atom.index for atom in residue.atoms()]
                expected_constraint_set.add(tuple(sorted([indices[2],
                                                          indices[0]])))
                expected_constraint_set.add(tuple(sorted([indices[2],
                                                          indices[1]])))   
                expected_constraint_set.add(tuple(sorted([indices[1],
                                                          indices[0]])))
    if constraints_value==HBonds:
        for bond in topology.bonds():
            if bond[0].element.symbol=='H' or bond[1].element.symbol=='H':
                expected_constraint_set.add(tuple(sorted([bond[0].index, 
                                                          bond[1].index])))
    if constraints_value in set([AllBonds, HAngles]):
        bonds=[b for b in topology.bonds()]
        for b in bonds:
            expected_constraint_set.add(tuple(sorted([b[0].index, 
                                                      b[1].index])))
    if constraints_value==HAngles:
        expected_constraint_set=expected_constraint_set.union(constraintsHAngles(bonds))

    # Check that the size of the expected constraint set is the same as the
    # number of actual constraints in the system.
    self.assertEqual(len(expected_constraint_set), system.getNumConstraints())
    
    # Check each constraint in the system to see that it appears in the 
    # expected constraint set
    constraints = [system.getConstraintParameters(i) for i 
                   in range(system.getNumConstraints())]
    for c in constraints:
        self.assertTrue((c[0],c[1]) in expected_constraint_set or 
                        (c[1],c[0]) in expected_constraint_set)

    # Build the list of bonds in the topology.
    tbonds = [tuple(sorted([bond[0].index, bond[1].index])) for bond in topology.bonds()]

    # Build set of harmonic bond forces in the system.
    harmonic_bond_set = set()
    forces = system.getForces()
    for f in forces:
        if isinstance(f, HarmonicBondForce):
            for i in range(f.getNumBonds()):
                bond = f.getBondParameters(i)
                harmonic_bond_set.add(tuple(sorted([bond[0],bond[1]])))

    # Build set of constraints in the system.
    constraints = [system.getConstraintParameters(i) for i 
                   in range(system.getNumConstraints())]
    constraint_set = set()
    for c in constraints:
        constraint_set.add(tuple(sorted([c[0], c[1]])))
        
    # Check that each bond in the topology is either in the harmonic bond set or in 
    # the constraint set, but not both.
    for bond in tbonds:
        self.assertNotEqual(bond in harmonic_bond_set, bond in constraint_set)

def constraintsHAngles(bonds):
    """Given a set of bonds from a Topology, return the set of angle constraints that 
    you expect to appear in the system.  An angle is constrained if the bond sequence 
    is H-O-X or H-X-H where X is any element.

    """
    expected_constraint_set = set()

    for bond in bonds:
        if bond[0].element.symbol=='H' and bond[1].element.symbol=='O':
            indexH = bond[0].index
            indexO = bond[1].index
            for bond2 in bonds:
                if bond2[0].index==indexO and bond2[1].index!=indexH:
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[1].index])))
                elif bond2[1].index==indexO and bond2[0].index!=indexH:
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[0].index])))
        elif bond[1].element.symbol=='H' and bond[0].element.symbol=='O':
            indexH = bond[1].index
            indexO = bond[0].index
            for bond2 in bonds:
                if bond2[0].index==indexO and bond2[1].index!=indexH:
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[1].index])))
                elif bond2[1].index==indexO and bond2[0].index!=indexH:
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[0].index])))
        elif bond[0].element.symbol=='H' and bond[1].element.symbol!='O':
            indexH = bond[0].index
            indexX = bond[1].index
            for bond2 in bonds:
                if bond2[0].index==indexX and (bond2[1].index!=indexH 
                                               and bond2[1].element.symbol=='H'):
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[1].index])))
                elif bond2[1].index==indexX and (bond2[0].index!=indexH 
                                                 and bond2[0].element.symbol=='H'):
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[0].index])))
        elif bond[1].element.symbol=='H' and bond[0].element.symbol!='O':
            indexH = bond[1].index
            indexX = bond[0].index
            for bond2 in bonds:
                if bond2[0].index==indexX and (bond2[1].index!=indexH 
                                               and bond2[1].element.symbol=='H'):
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[1].index])))
                elif bond2[1].index==indexX and (bond2[0].index!=indexH 
                                                 and bond2[0].element.symbol=='H'):
                    expected_constraint_set.add(tuple(sorted([indexH, bond2[0].index])))

    return expected_constraint_set
