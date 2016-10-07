"""
Convert RCSB components (monomers) into residues.xml file.

Requires 'components-pub-xml.tar.gz' from http://ligand-expo.rcsb.org/dictionaries/components-pub-xml.tar.gz
"""
import os
import lxml.etree as etree
from gzip import GzipFile
import tarfile

addTerminal = ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','CYX','PRO','PHE',
                   'TYR','TRP','HIS','ASP','ASN','GLU','GLN','MET','LYS',
                   'ARG','ORN','AIB','PCA','FOR','UNK']
addProtBackbone = ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','CYX','PRO','PHE',
                   'TYR','TRP','HIS','ASP','ASN','GLU','GLN','MET','LYS',
                   'ARG','ORN','AIB','PCA','FOR','NME','NH2','UNK']
addNucBackbone = ['A','G','C','U','DA','DG','DC','DT']

# All amino acids have these atoms
amino_acid_atoms = set(['N', 'CA', 'C', 'O'])

# All nucleic acids have these atoms
nucleic_acid_atoms = set(['OP3', 'P'])

def parse_monomer(infile, outfile):
    parser = etree.XMLParser(remove_blank_text=True) # For pretty print on write
    tree = etree.parse(infile, parser)
    root = tree.getroot()

    # Write residue name
    resname = root.attrib['datablockName']
    outfile.write(' <Residue name="%s">\n' % resname)

    atoms = set()
    bonds = set()
    for category in root:
        # Extract atoms
        if category.tag == '{http://pdbml.pdb.org/schema/pdbx-v40.xsd}chem_comp_atomCategory':
            for atom_node in category:
                atom = atom_node.attrib['atom_id']
                atoms.add(atom)
        # Extract bonds
        elif category.tag == '{http://pdbml.pdb.org/schema/pdbx-v40.xsd}chem_comp_bondCategory':
            for bond_node in category:
                atom1 = bond_node.attrib['atom_id_1']
                atom2 = bond_node.attrib['atom_id_2']
                bonds.add((min(atom1, atom2), max(atom1, atom2)))

    # Add protein backbone and terminal bonds
    if amino_acid_atoms.issubset(atoms):
        bonds.add(("-C", "N"))
        bonds.add(("H3", "N"))
    # Add nucleic acid bonds
    elif nucleic_acid_atoms.issubset(atoms):
        bonds.add(("-OP3", "P"))

    # Write bonds
    for bond in sorted(bonds):
        outfile.write('  <Bond from="%s" to="%s"/>\n' % bond)
    outfile.write(' </Residue>\n')

#
# Generate residues file by processing components-pub-xml
#

# Open output file.
outfile = GzipFile('residues.xml.gz', 'w:gz')

# Write header
outfile.write("<Residues>\n")

# Parse all monomers
filename = 'components-pub-xml.tar.gz'
tarball = tarfile.open(filename, mode='r:gz')
for member in tarball:
    if member.isfile():
        print(member)
        # Extract file
        file = tarball.extractfile(member)
        parse_monomer(file, outfile)

# Write footer
outfile.write("</Residues>\n")
