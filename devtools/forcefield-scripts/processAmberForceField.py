import sys
import math
import simtk.openmm.app.element as element
import simtk.unit as unit

elements = {}
for elem in element.Element._elements_by_symbol.values():
    num = elem.atomic_number
    if num not in elements or elem.mass < elements[num].mass:
        elements[num] = elem

OTHER = 0
ATOMS = 1
CONNECT = 2
CONNECTIVITY = 3
RESIDUECONNECT = 4
section = OTHER

residueAtoms = {}
residueBonds = {}
residueConnections = {}

types = []
masses = {}
resAtomTypes = {}
vdwEquivalents = {}
vdw = {}
charge = {}
bonds = []
angles = []
torsions = []
impropers = []
charge14scale = 1.0/1.2
epsilon14scale = 0.5

skipResidues = ['CIO', 'IB'] # "Generic" ions defined by Amber, which are identical to other real ions
skipClasses = ['OW', 'HW'] # Skip water atoms, since we define these in separate files

def addAtom(residue, atomName, atomClass, element, charge):
    if residue is None:
        return
    residueAtoms[residue].append([atomName, len(types)])
    types.append((atomClass, element, charge))

def addBond(residue, atom1, atom2):
    if residue is None:
        return
    residueBonds[residue].append((atom1, atom2))

def addExternalBond(residue, atom):
    if residue is None:
        return
    if atom != -1:
        residueConnections[residue] += [atom]

# Load input files.

for inputfile in sys.argv[1:]:
    if inputfile.endswith('.lib') or inputfile.endswith('.off'):
        # Read a library file
        for line in open(inputfile):
            if line.startswith('!entry'):
                fields = line.split('.')
                residue = fields[1]
                if residue in skipResidues:
                    residue = None
                    continue
                key = fields[3].split()[0]
                if key == 'atoms':
                    section = ATOMS
                    residueAtoms[residue] = []
                    residueBonds[residue] = []
                    residueConnections[residue] = []
                elif key == 'connect':
                    section = CONNECT
                elif key == 'connectivity':
                    section = CONNECTIVITY
                elif key == 'residueconnect':
                    section = RESIDUECONNECT
                else:
                    section = OTHER
            elif section == ATOMS:
                fields = line.split()
                atomName = fields[0][1:-1]
                atomClass = fields[1][1:-1]
                if fields[6] == '-1':
                    # Workaround for bug in some Amber files.
                    if atomClass[0] == 'C':
                        elem = elements[6]
                    elif atomClass[0] == 'H':
                        elem = elements[1]
                    else:
                        raise ValueError('Illegal atomic number: '+line)
                else:
                    elem = elements[int(fields[6])]
                charge = float(fields[7])
                addAtom(residue, atomName, atomClass, elem, charge)
            elif section == CONNECT:
                addExternalBond(residue, int(line)-1)
            elif section == CONNECTIVITY:
                fields = line.split()
                addBond(residue, int(fields[0])-1, int(fields[1])-1)
            elif section == RESIDUECONNECT:
                # Some Amber files have errors in them, incorrectly listing atoms that should not be
                # connected in the first two positions.  We therefore rely on the "connect" section for
                # those, using this block only for other external connections.
                for atom in [int(x)-1 for x in line.split()[2:]]:
                    addExternalBond(residue, atom)

    elif inputfile.endswith('.in'):
        lines = open(inputfile).read().split('\n')
        i = 2
        while True:
            if lines[i].strip() == 'STOP':
                break
            i += 1
            residue = lines[i].strip()
            # Hack to get unique residue names, since different files use the same names.
            if inputfile.endswith('nt.in'):
                residue = 'N'+residue
            if inputfile.endswith('ct.in'):
                residue = 'C'+residue
            residueAtoms[residue] = []
            residueBonds[residue] = []
            residueConnections[residue] = []
            atoms = []
            mainchain = []
            i += 7
            while len(lines[i].rstrip()) > 0:
                fields = lines[i].split()
                atoms.append((fields[1], fields[2], float(fields[10]))) # (name, type, charge)
                if fields[3] == 'M':
                    mainchain.append(len(atoms)-1)
                bondedTo = int(fields[4])-4
                if bondedTo >= 0:
                    addBond(residue, len(atoms)-1, bondedTo)
                i += 1
            while True:
                if lines[i].strip() == 'LOOP':
                    i += 1
                    while len(lines[i].rstrip()) > 0:
                        fields = lines[i].strip().split()
                        bondFrom = [j for j in range(len(atoms)) if fields[0] == atoms[j][0]][0]
                        bondTo = [j for j in range(len(atoms)) if fields[1] == atoms[j][0]][0]
                        addBond(residue, bondFrom, bondTo)
                        i += 1
                if lines[i].strip() != 'DONE':
                    i += 1
                    continue
                for atom in atoms:
                    try:
                        el = element.Element.getBySymbol(atom[1][0])
                    except:
                        el = None
                    addAtom(residue, atom[0], atom[1], el, atom[2])
                # Hack for figuring out external bonds.  We'll have to fix up disulfide bonds and capping groups manually.
                if len(mainchain) > 0 and not inputfile.endswith('nt.in'):
                    addExternalBond(residue, mainchain[0])
                if len(mainchain) > 1 and not inputfile.endswith('ct.in'):
                    addExternalBond(residue, mainchain[-1])
                i += 1
                break
            
    
    elif inputfile.endswith('.dat'):
        # Read a force field file.
        block = 0
        continueTorsion = False
        for line in open(inputfile):
            line = line.strip()
            if block == 0:     # Title
                block += 1
            elif block == 1:   # Mass
                fields = line.split()
                if len(fields) == 0:
                    block += 1
                else:
                    masses[fields[0]] = float(fields[1])
            elif block == 2:   # Hydrophilic atoms
                block += 1
            elif block == 3:   # Bonds
                if len(line) == 0:
                    block += 1
                elif '-' in line:
                    fields = line[5:].split()
                    bonds.append((line[:2].strip(), line[3:5].strip(), fields[0], fields[1]))
            elif block == 4:   # Angles
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[8:].split()
                    angles.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), fields[0], fields[1]))
            elif block == 5:   # Torsions
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[11:].split()
                    periodicity = int(float(fields[3]))
                    if continueTorsion:
                        torsions[-1] += [float(fields[1])/float(fields[0]), fields[2], abs(periodicity)]
                    else:
                        torsions.append([line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), float(fields[1])/float(fields[0]), fields[2], abs(periodicity)])
                    continueTorsion = (periodicity < 0)
            elif block == 6:   # Improper torsions
                if len(line) == 0:
                    block += 1
                else:
                    fields = line[11:].split()
                    impropers.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), fields[0], fields[1], fields[2]))
            elif block == 7:   # 10-12 hbond potential
                if len(line) == 0:
                    block += 1
            elif block == 8:   # VDW equivalents
                if len(line) == 0:
                    block += 1
                else:
                    fields = line.split()
                    for atom in fields[1:]:
                        vdwEquivalents[atom] = fields[0]
            elif block == 9:   # VDW type
                block += 1
                vdwType = line.split()[1]
                if vdwType not in ['RE', 'AC']:
                    raise ValueError('Nonbonded type (KINDNB) must be RE or AC') 
            elif block == 10:   # VDW parameters
                if len(line) == 0:
                    block += 1
                else:
                    fields = line.split()
                    vdw[fields[0]] = (fields[1], fields[2])
    
    else:
        # Assume it's a frcmod file.
        block = ''
        continueTorsion = False
        first = True
        for line in open(inputfile):
            line = line.strip()
            if len(line) == 0 or first:
                block = None
                first = False
            elif block is None:
                block = line
            elif block.startswith('MASS'):
                fields = line.split()
                masses[fields[0]] = float(fields[1])
            elif block.startswith('BOND'):
                fields = line[5:].split()
                bonds.append((line[:2].strip(), line[3:5].strip(), fields[0], fields[1]))
            elif block.startswith('ANGL'):
                fields = line[8:].split()
                angles.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), fields[0], fields[1]))
            elif block.startswith('DIHE'):
                fields = line[11:].split()
                periodicity = int(float(fields[3]))
                if continueTorsion:
                    torsions[-1] += [float(fields[1])/float(fields[0]), fields[2], abs(periodicity)]
                else:
                    torsions.append([line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), float(fields[1])/float(fields[0]), fields[2], abs(periodicity)])
                continueTorsion = (periodicity < 0)
            elif block.startswith('IMPR'):
                fields = line[11:].split()
                impropers.append((line[:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip(), fields[0], fields[1], fields[2]))
            elif block.startswith('NONB'):
                fields = line.split()
                vdw[fields[0]] = (fields[1], fields[2])

# Reduce the list of atom types.  If multiple hydrogens are bound to the same heavy atom,
# they should all use the same type.

removeType = [False]*len(types)
for res in residueAtoms:
    if res not in residueBonds:
        continue
    atomBonds = [[] for atom in residueAtoms[res]]
    for bond in residueBonds[res]:
        atomBonds[bond[0]].append(bond[1])
        atomBonds[bond[1]].append(bond[0])
    for index, atom in enumerate(residueAtoms[res]):
        hydrogens = [x for x in atomBonds[index] if types[residueAtoms[res][x][1]][1] == element.hydrogen]
        for h in hydrogens[1:]:
            removeType[residueAtoms[res][h][1]] = True
            residueAtoms[res][h][1] = residueAtoms[res][hydrogens[0]][1]
newTypes = []
replaceWithType = [0]*len(types)
for i in range(len(types)):
    if not removeType[i]:
        newTypes.append(types[i])
    replaceWithType[i] = len(newTypes)-1
types = newTypes
for res in residueAtoms:
    for atom in residueAtoms[res]:
        atom[1] = replaceWithType[atom[1]]

# Create the XML output.

def fix(atomClass):
    if atomClass == 'X':
        return ''
    return atomClass

print "<ForceField>"
print " <AtomTypes>"
for index, type in enumerate(types):
    if type[1] is None:
        el = ""
        mass = 0
    else:
        el = type[1].symbol
        mass = type[1].mass.value_in_unit(unit.amu)
    print """  <Type name="%d" class="%s" element="%s" mass="%s"/>""" % (index, type[0], el, mass)
print " </AtomTypes>"
print " <Residues>"
for res in sorted(residueAtoms):
    print """  <Residue name="%s">""" % res
    for atom in residueAtoms[res]:
        print "   <Atom name=\"%s\" type=\"%d\"/>" % tuple(atom)
    if res in residueBonds:
        for bond in residueBonds[res]:
            print """   <Bond from="%d" to="%d"/>""" % bond
    if res in residueConnections:
        for bond in residueConnections[res]:
            print """   <ExternalBond from="%d"/>""" % bond
    print "  </Residue>"
print " </Residues>"
print " <HarmonicBondForce>"
processed = set()
for bond in bonds:
    signature = (bond[0], bond[1])
    if signature in processed:
        continue
    if any([c in skipClasses for c in signature]):
        continue
    processed.add(signature)
    length = float(bond[3])*0.1
    k = float(bond[2])*2*100*4.184
    print """  <Bond class1="%s" class2="%s" length="%s" k="%s"/>""" % (bond[0], bond[1], str(length), str(k))
print " </HarmonicBondForce>"
print " <HarmonicAngleForce>"
processed = set()
for angle in angles:
    signature = (angle[0], angle[1], angle[2])
    if signature in processed:
        continue
    if any([c in skipClasses for c in signature]):
        continue
    processed.add(signature)
    theta = float(angle[4])*math.pi/180.0
    k = float(angle[3])*2*4.184
    print """  <Angle class1="%s" class2="%s" class3="%s" angle="%s" k="%s"/>""" % (angle[0], angle[1], angle[2], str(theta), str(k))
print " </HarmonicAngleForce>"
print " <PeriodicTorsionForce>"
processed = set()
for tor in reversed(torsions):
    signature = (fix(tor[0]), fix(tor[1]), fix(tor[2]), fix(tor[3]))
    if signature in processed:
        continue
    if any([c in skipClasses for c in signature]):
        continue
    processed.add(signature)
    tag = "  <Proper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\"" % signature
    i = 4
    while i < len(tor):
        index = i/3
        periodicity = int(float(tor[i+2]))
        phase = float(tor[i+1])*math.pi/180.0
        k = tor[i]*4.184
        tag += " periodicity%d=\"%d\" phase%d=\"%s\" k%d=\"%s\"" % (index, periodicity, index, str(phase), index, str(k))
        i += 3
    tag += "/>"
    print tag
processed = set()
for tor in reversed(impropers):
    signature = (fix(tor[2]), fix(tor[0]), fix(tor[1]), fix(tor[3]))
    if signature in processed:
        continue
    if any([c in skipClasses for c in signature]):
        continue
    processed.add(signature)
    tag = "  <Improper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\"" % signature
    i = 4
    while i < len(tor):
        index = i/3
        periodicity = int(float(tor[i+2]))
        phase = float(tor[i+1])*math.pi/180.0
        k = float(tor[i])*4.184
        tag += " periodicity%d=\"%d\" phase%d=\"%s\" k%d=\"%s\"" % (index, periodicity, index, str(phase), index, str(k))
        i += 3
    tag += "/>"
    print tag
print " </PeriodicTorsionForce>"
print """ <NonbondedForce coulomb14scale="%g" lj14scale="%s">""" % (charge14scale, epsilon14scale)
sigmaScale = 0.1*2.0/(2.0**(1.0/6.0))
for index, type in enumerate(types):
    atomClass = type[0]
    q = type[2]
    if atomClass in vdwEquivalents:
        atomClass = vdwEquivalents[atomClass]
    if atomClass in vdw:
        params = [float(x) for x in vdw[atomClass]]
        if vdwType == 'RE':
            sigma = params[0]*sigmaScale
            epsilon = params[1]*4.184
        else:
            sigma = (params[0]/params[1])**(1.0/6.0)
            epsilon = 4.184*params[1]*params[1]/(4*params[0])
    else:
        sigma = 1
        epsilon = 0
    if sigma == 0 or epsilon == 0:
        sigma, epsilon = 1, 0
    if q != 0 or epsilon != 0:
        print """  <Atom type="%d" charge="%s" sigma="%s" epsilon="%s"/>""" % (index, q, sigma, epsilon)
print " </NonbondedForce>"
print "</ForceField>"

