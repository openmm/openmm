import os

addTerminal = ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','CYX','PRO','PHE',
                   'TYR','TRP','HIS','ASP','ASN','GLU','GLN','MET','LYS',
                   'ARG','ORN','AIB','PCA','FOR','UNK']
addProtBackbone = ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','CYX','PRO','PHE',
                   'TYR','TRP','HIS','ASP','ASN','GLU','GLN','MET','LYS',
                   'ARG','ORN','AIB','PCA','FOR','NME','NH2','UNK']
addNucBackbone = ['A','G','C','U','DA','DG','DC','DT']
print "<Residues>"
for file in os.listdir('monomers'):
    try:
        f = open('monomers/'+file)
    except IOError:
        continue
    print """ <Residue name="%s">""" % file
    bonds = set()
    for line in f:
        fields = line.split()
        if fields[0] == 'CONECT':
            a1 = fields[1]
            for a2 in fields[3:]:
                bonds.add((min(a1, a2), max(a1, a2)))
    if file in addProtBackbone:
        bonds.add(("-C", "N"))
    elif file in addNucBackbone:
        bonds.add(("-OP3", "P"))
    if file in addTerminal:
        bonds.add(("H3", "N"))        
    for bond in sorted(bonds):
        print """  <Bond from="%s" to="%s"/>""" % bond
    print " </Residue>"
print "</Residues>"

