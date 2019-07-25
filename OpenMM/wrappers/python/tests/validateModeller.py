
def validate_preserved(self, topology_before, topology_after, chain_dict, residue_dict, atom_dict):
    """ Validate that the residues and atoms are correctly preserved in the given topologies. """    

    # validate that residues are preserved
    residues_before = [residue for residue in topology_before.residues()]
    residues_after = [residue for residue in topology_after.residues()]
    for residue_index in residue_dict:
        # validate that before and after residues have the same name
        self.assertTrue(residues_before[residue_index].name == 
                        residues_after[residue_dict[residue_index]].name)
        # validate that before and after residues come from the same chain
        self.assertTrue(chain_dict[residues_before[residue_index].chain.index] ==
                        residues_after[residue_dict[residue_index]].chain.index)

    # validate that atoms are preserved
    atoms_before = [atom for atom in topology_before.atoms()]
    atoms_after = [atom for atom in topology_after.atoms()]
    for atom_index in atom_dict:
        # validate that before and after atoms have the same name
        self.assertTrue(atoms_before[atom_index].name == 
                        atoms_after[atom_dict[atom_index]].name)
        # validate that before and after atoms are the same element
        self.assertTrue(atoms_before[atom_index].element ==
                        atoms_after[atom_dict[atom_index]].element)
        # validate that before and after atoms come from the same residue
        self.assertTrue(residue_dict[atoms_before[atom_index].residue.index] ==
                        atoms_after[atom_dict[atom_index]].residue.index)
# end of validate_preserved()

def validate_deltas(self, topology_before, topology_after, chain_delta, residue_delta, atoms_delta):
    """ Validate that the number of chains, residues, and atoms have changes the expected amounts. """

    # validate that the correct number of chains are deleted
    chains_before = [chain for chain in topology_before.chains()]
    chains_after = [chain for chain in topology_after.chains()]
    self.assertTrue(len(chains_after) == len(chains_before) + chain_delta)

    # validate that the correct number of residues are deleted
    residues_before = [residue for residue in topology_before.residues()]
    residues_after = [residue for residue in topology_after.residues()]
    self.assertTrue(len(residues_after) == len(residues_before) + residue_delta)

    # validate that the correct number of atoms are deleted
    atoms_before = [atom for atom in topology_before.atoms()]
    atoms_after = [atom for atom in topology_after.atoms()]
    self.assertTrue(len(atoms_after) == len(atoms_before) + atoms_delta)
# end of validate_deltas()

def validate_equivalence(self, topology_before, topology_after):
    """ Validate that two topologies are equivalent. """

    # First, check that they have the same number of each type of atom.

    # construct dictionaries with the counts of each atom by symbol
    atom_count_before = {}
    for atom in topology_before.atoms():
        if atom.element.symbol not in atom_count_before:
            atom_count_before[atom.element.symbol] = 1
        else:
            atom_count_before[atom.element.symbol] += 1

    atom_count_after = {}
    for atom in topology_after.atoms():
        if atom.element.symbol not in atom_count_after:
            atom_count_after[atom.element.symbol] = 1
        else:
            atom_count_after[atom.element.symbol] += 1
    
    # make sure the dictionaries are equivalent
    self.assertTrue(atom_count_before==atom_count_after)
    
    # Next, we want to confirm that every hydrogen is in exactly one bond
    # and that heavy atoms are bonded to the same number of hydrogens.
    # To achieve this, we build two dictionaries for each topology; one for 
    # hydrogens and one for heavy atoms.  The keys will be the atom indices.

    def build_bond_dicts(topology):
        hydrogen_bonds = {}
        heavy_atom_bonds = {}    
        for bond in topology.bonds():
            if bond[0].element.symbol=='H':
                if bond[0].index not in hydrogen_bonds:
                    hydrogen_bonds[bond[0].index] = 1
                else:
                    hydrogen_bonds[bond[0].index] += 1
                if bond[1].index not in heavy_atom_bonds:
                    heavy_atom_bonds[bond[1].index] = 1
                else:
                    heavy_atom_bonds[bond[1].index] += 1
            elif bond[1].element.symbol=='H':
                if bond[1].index not in hydrogen_bonds:
                    hydrogen_bonds[bond[1].index] = 1
                else:
                    hydrogen_bonds[bond[1].index] += 1
                if bond[0].index not in heavy_atom_bonds:
                    heavy_atom_bonds[bond[0].index] = 1
                else:
                    heavy_atom_bonds[bond[0].index] += 1    
        
        return (hydrogen_bonds, heavy_atom_bonds)

    (hydrogen_bonds_before, heavy_atom_bonds_before) = build_bond_dicts(topology_before)  
    (hydrogen_bonds_after, heavy_atom_bonds_after) = build_bond_dicts(topology_after)      

    atom_list_before = []
    for atom in topology_before.atoms():
        if atom.element.symbol != 'H':
            atom_list_before.append((atom.index, atom.name))
    atom_list_after = []
    for atom in topology_after.atoms():
        if atom.element.symbol != 'H':
            atom_list_after.append((atom.index, atom.name))

    correspondence_dict = {}
    for i in range(len(atom_list_before)):
        correspondence_dict[atom_list_before[i][0]] = atom_list_after[i][0]
        if atom_list_before[i][1] != atom_list_after[i][1]:
            # shouldn't get here!
            self.assertTrue(0==1)

    # also check that there is exactly one bond per hydrogen after adding hydrogens
    for key in hydrogen_bonds_after:
        self.assertTrue(hydrogen_bonds_after[key]==1)

    # check that each heavy atom has the same number of hydrogen bonds
    for key in heavy_atom_bonds_before:
        self.assertTrue(heavy_atom_bonds_before[key]==heavy_atom_bonds_after[correspondence_dict[key]])

# end of validate_equivalence
