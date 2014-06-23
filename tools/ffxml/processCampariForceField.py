import sys
import math
import shlex
import xml.etree.ElementTree as etree
import simtk.pyopenmm.element as element
import simtk.unit as unit
from simtk.openmm.app import PDBFile

elements = {}
for elem in element.Element.elements_by_symbol.values():
    num = elem.atomic_number
    if num not in elements or elem.mass < elements[num].mass:
        elements[num] = elem
PDBFile._loadNameReplacementTables()

residueNames = {
    "Acetyl N-Terminus":"ACE",
    "Alanine":"ALA",
    "Arginine(+)":"ARG",
    "Asparagine":"ASN",
    "Aspartic Acid(-)":"ASP",
    "C-Terminal Glycine(-)":"GLY",
    "C-Terminal Proline(-)":"PRO",
    "C-Terminal Residue(-)":"Protein",
    "Cysteine (-SH)":"CYS",
    "Cystine (-SS-)":"CYX",
    "Glutamic Acid(-)":"GLU",
    "Glutamine":"GLN",
    "Glycine":"GLY",
    #"Histidine (+)":"HIP",
    "Histidine (HD)":"HID",
    "Histidine (HE)":"HIE",
    "Isoleucine":"ILE",
    "Leucine":"LEU",
    "Lysine(+)":"LYS",
    "Methionine":"MET",
    "N-MeAmide C-Terminus":"NME",
    "N-Terminal Glycine(+)":"GLY",
    "N-Terminal Proline(+)":"PRO",
    "N-Terminal Residue(+)":"Protein",
    "Phenylalanine":"PHE",
    "Proline":"PRO",
    "Serine":"SER",
    "Threonine":"THR",
    "Tryptophan":"TRP",
    "Tyrosine":"TYR",
    "Valine":"VAL",
    "Bromide Ion":"Br-",
    "Caesium Ion":"Cs+",
    "Chloride Ion":"Cl-",
    "Iodide Ion":"I-",
    "Potassium Ion":"K+",
    "Sodium Ion":"Na+"
}

# Load the force field file.

biotypes = []
atomTypes = {}
sigma = {}
epsilon = {}
charge = {}
radius = {}
residueBiotypes = {}
nResidueBiotypes = {}
cResidueBiotypes = {}
bondTypes = {}
angleTypes = {}
torsionTypes = {}
bonds = []
angles = []
torsions = []
impropers = []
maxLJType = 0
for res in residueNames.itervalues():
    residueBiotypes[res] = {}
    nResidueBiotypes[res] = {}
    cResidueBiotypes[res] = {}

for line in open(sys.argv[1]):
    try:
        fields = shlex.split(line)
        if len(fields) == 0:
            continue
        if fields[0] == 'biotype':
            if '0' not in fields[4:6] and fields[3] in residueNames:
                biotypes.append(fields[1:])
                res = residueNames[fields[3]]
                if res in ['HID', 'HIE', 'HIP']:
                    standardName = 'HIS'
                elif res == 'CYX':
                    standardName = 'CYS'
                else:
                    standardName = res
                atomName = fields[2]
                if standardName in PDBFile._atomNameReplacements and atomName in PDBFile._atomNameReplacements[standardName]:
                    atomName = PDBFile._atomNameReplacements[standardName][atomName]
                if fields[3].startswith('N-Terminal'):
                    nResidueBiotypes[res][atomName] = fields[1]
                elif fields[3].startswith('C-Terminal'):
                    cResidueBiotypes[res][atomName] = fields[1]
                else:
                    residueBiotypes[res][atomName] = fields[1]
        elif fields[0] == 'atom':
            atomTypes[fields[1]] = (elements[int(fields[4])], fields[5])
        elif fields[0] == 'contact':
            type1 = int(fields[1])-1
            type2 = int(fields[2])-1
            if type1 not in sigma:
                sigma[type1] = {}
            if type2 not in sigma:
                sigma[type2] = {}
            value = float(fields[3])/10
            sigma[type1][type2] = value
            sigma[type2][type1] = value
            maxLJType = max(maxLJType, type1, type2)
        elif fields[0] == 'interact':
            if fields[1] == fields[2]:
                epsilon[fields[1]] = float(fields[3])*4.184
        elif fields[0] == 'radius':
            radius[fields[1]] = float(fields[2])
        elif fields[0] == 'charge':
            charge[fields[1]] = float(fields[3])
        elif fields[0] == 'bond':
            if fields[2] == '1':
                bondTypes[fields[1]] = (float(fields[3])*2*418.4, float(fields[4])/10)
        elif fields[0] == 'angle':
            if fields[2] == '1':
                angleTypes[fields[1]] = (float(fields[3])*2*4.184, float(fields[4])*math.pi/180.0)
        elif fields[0] == 'torsion':
            if fields[2] == '1':
                params = [float(fields[i+3])*4.184*(-1)**i for i in range(6)]
                torsionTypes[fields[1]] = [0 if x == -0 else x for x in params]
        elif fields[0] == 'bonded_type_bond':
            bonds.append(fields[1:])
        elif fields[0] == 'bonded_type_angle':
            angles.append(fields[1:])
        elif fields[0] == 'bonded_type_torsion':
            torsions.append(fields[1:])
        elif fields[0] == 'bonded_type_imptors':
            impropers.append(fields[1:])
    except:
        pass

for type in epsilon:
    if type not in radius:
        intType = int(type)-1
        radius[type] = sigma[intType][intType]/2
for res in residueNames.itervalues():
    for name in nResidueBiotypes['Protein']:
        if name not in nResidueBiotypes[res]:
            nResidueBiotypes[res][name] = nResidueBiotypes['Protein'][name]
    for name in cResidueBiotypes['Protein']:
        if name not in cResidueBiotypes[res]:
            cResidueBiotypes[res][name] = cResidueBiotypes['Protein'][name]
    for name in residueBiotypes[res]:
        if name not in nResidueBiotypes[res]:
            nResidueBiotypes[res][name] = residueBiotypes[res][name]
        if name not in cResidueBiotypes[res]:
            cResidueBiotypes[res][name] = residueBiotypes[res][name]

# Load and process parameters the were extracted directly from CAMPARI.

volumeParams = {}
for line in open('absinthVolumeParams.txt'):
    fields = shlex.split(line)
    type = fields[1]
    maxsav = float(fields[2])
    volReduction = float(fields[3])
    if type not in volumeParams:
        volumeParams[type] = (1, maxsav, volReduction)
    else:
        values = volumeParams[type]
        volumeParams[type] = (values[0]+1, values[1]+maxsav, values[2]+volReduction)
solvationParams = {}
groups = {}
for line in open('absinthSolvationGroups.txt'):
    fields = shlex.split(line)
    residue = int(fields[0])
    type = fields[3]
    group = int(fields[1])
    fos = float(fields[5])
    solvationParams[type] = fos
    groupKey = (residue, group)
    if groupKey in groups:
        groups[groupKey].add(type)
    else:
        groups[groupKey] = set([type])
uniqueGroups = set(tuple(val) for val in groups.itervalues())
solvationGroups = {}
for type in biotypes:
    if type[1] in ('C', 'O'):
        solvationGroups[type[0]] = 0
    elif type[1] in ('N', 'HN'):
        solvationGroups[type[0]] = 1
nextSolvGroupId = 2
for group in uniqueGroups:
    id = None
    for type in group:
        if type in solvationGroups:
            id = solvationGroups[type]
    if id is None:
        id = nextSolvGroupId
        nextSolvGroupId += 1
    for type in group:
        solvationGroups[type] = id
groups = {}
for line in open('absinthChargeGroups.txt'):
    fields = shlex.split(line)
    residue = int(fields[0])
    type = fields[3]
    group = int(fields[1])
    groupKey = (residue, group)
    if groupKey in groups:
        groups[groupKey].add(type)
    else:
        groups[groupKey] = set([type])
uniqueGroups = set(tuple(val) for val in groups.itervalues())
chargeGroups = {}
nextChargeGroupId = 0
for group in uniqueGroups:
    id = None
    for type in group:
        if type in chargeGroups:
            id = chargeGroups[type]
    if id is None:
        id = nextChargeGroupId
        nextChargeGroupId += 1
    for type in group:
        chargeGroups[type] = id

# Generate the XML file.

print '<ForceField>'
print ' <AtomTypes>'
for type in biotypes:
    print '  <Type name="%s" class="%s" element="%s" mass="%s"/>' % (type[0], type[5], atomTypes[type[3]][0].symbol, atomTypes[type[3]][1])
print ' </AtomTypes>'

def findType(name, residue, types):
    while len(name) > 0:
        if name in types[residue]:
            return types[residue][name]
        if name in types['Protein']:
            return types['Protein'][name]
        if name[-1].isdigit() or name[-1] == '+' or name[-1] == '-':
            name = name[:-1]
        else:
            return None
    return None

typeEquivalents = {}
tree = etree.parse(sys.argv[2])
print ' <Residues>'
for residue in tree.getroot().find('Residues').findall('Residue'):
    resName = residue.attrib['name']
    if len(resName) == 4 and resName[0] == 'N':
        types = nResidueBiotypes
        shortName = resName[1:]
    elif len(resName) == 4 and resName[0] == 'C':
        types = cResidueBiotypes
        shortName = resName[1:]
    else:
        types = residueBiotypes
        shortName = resName
    if shortName not in types:
        continue
    print '  <Residue name="%s">' % resName
    for atom in residue.findall('Atom'):
        atomName = atom.attrib['name'].upper()
        if shortName in PDBFile._atomNameReplacements and atomName in PDBFile._atomNameReplacements[shortName]:
            atomName = PDBFile._atomNameReplacements[shortName][atomName]
        type = findType(atomName, shortName, types)
        if type is None and atom.attrib['type'] in typeEquivalents:
            type = typeEquivalents[atom.attrib['type']]
        if type is None:
            print atomName, atom.attrib['name'], resName, 'None'
        else:
            typeEquivalents[atom.attrib['type']] = type
            print '   <Atom name="%s" type="%s"/>' % (atomName, type)
    for bond in residue.findall('Bond'):
        print '   <Bond from="%s" to="%s"/>' % (bond.attrib['from'], bond.attrib['to'])
    for bond in residue.findall('ExternalBond'):
        print '   <ExternalBond from="%s"/>' % bond.attrib['from']
    print '  </Residue>'
print ' </Residues>'
print ' <HarmonicBondForce>'
for bond in bonds:
    if bond[2] in bondTypes:
        type = bondTypes[bond[2]]
        print '  <Bond class1="%s" class2="%s" length="%g" k="%g"/>' % (bond[0], bond[1], type[1], type[0])
print ' </HarmonicBondForce>'
print ' <HarmonicAngleForce>'
for angle in angles:
    if angle[3] in angleTypes:
        type = angleTypes[angle[3]]
        print '  <Angle class1="%s" class2="%s" class3="%s" angle="%g" k="%g"/>' % (angle[0], angle[1], angle[2], type[1], type[0])
print ' </HarmonicAngleForce>'
print ' <RBTorsionForce>'
for torsion in torsions:
    if torsion[4] in torsionTypes:
        type = torsionTypes[torsion[4]]
        print '  <Proper class1="%s" class2="%s" class3="%s" class4="%s" c0="%g" c1="%g" c2="%g" c3="%g" c4="%g" c5="%g"/>' % (torsion[0], torsion[1], torsion[2], torsion[3], type[0], type[1], type[2], type[3], type[4], type[5])
for torsion in impropers:
    if torsion[4] in torsionTypes:
        type = torsionTypes[torsion[4]]
        print '  <Improper class1="%s" class2="%s" class3="%s" class4="%s" c0="%g" c1="%g" c2="%g" c3="%g" c4="%g" c5="%g"/>' % (torsion[0], torsion[1], torsion[2], torsion[3], type[0], type[1], type[2], type[3], type[4], type[5])
print ' </RBTorsionForce>'
numLJTypes = maxLJType+1
print """ <CustomGBForce>
  <GlobalParameter name="numChargeGroups" defaultValue="0"/>
  <PerParticleParameter name="charge"/>
  <PerParticleParameter name="ljType"/>
  <PerParticleParameter name="epsilon"/>
  <PerParticleParameter name="radius"/>
  <PerParticleParameter name="chargeGroup"/>
  <PerParticleParameter name="solvGroup"/>
  <PerParticleParameter name="maxSav"/>
  <PerParticleParameter name="fos"/>
  <PerParticleParameter name="volReduction"/>
  <ComputedValue name="eta1" type="ParticlePairNoExclusions">
   gamma*volReduction2*(4*3.14159265/3)*radius2^3;
   gamma=anyOverlap*(outerOverlap*((dmax-r)/d2) + (1-outerOverlap)*(innerOverlap*(r/(radius1+radius2)) + (1-innerOverlap)));
   anyOverlap=step(dmax-r); outerOverlap=step(r-dmax+d2); innerOverlap=step(radius1+radius2-r);
   d2=2*radius2; dmax=radius1+radius2+rw; rw=0.5
  </ComputedValue>
  <ComputedValue name="eta" type="SingleParticle">
   1-eta1/vmax;
   vmax=(4*3.14159265/3)*((rw+radius)^3-radius^3); rw=0.5
  </ComputedValue>
  <ComputedValue name="screening" type="SingleParticle">
   1-a*v;
   v=max(0, min(1, d3+d2/(1+exp((d1-eta)/tau))));
   d3=1-d2/(1+exp((d1-maxSav)/tau));
   d2=1/(1/(1+exp((d1-maxSav)/tau)) - 1/(1+exp((d1-minSav)/tau)));
   d1=chi*maxSav+(1-chi)*minSav;
   minSav=0.2595; tau=0.5; chi=0.9; a=1-1/sqrt(78.2);
  </ComputedValue>
  <EnergyTerm type="SingleParticle">
   fos*v;
   v=max(0, min(1, d3+d2/(1+exp((d1-eta)/tau))));
   d3=1-d2/(1+exp((d1-maxSav)/tau));
   d2=1/(1/(1+exp((d1-maxSav)/tau)) - 1/(1+exp((d1-minSav)/tau)));
   d1=chi*maxSav+(1-chi)*minSav;
   minSav=0.2595; tau=0.25; chi=0.1
  </EnergyTerm>
  <EnergyTerm type="ParticlePair">
   4*eps*((sig/r)^12-(sig/r)^6); sig=sigma(ljType1*%d+ljType2); eps=sqrt(epsilon1*epsilon2)
  </EnergyTerm>
  <EnergyTerm type="ParticlePairNoExclusions">
   include(chargeGroup1*numChargeGroups+chargeGroup2)*screening1*screening2*138.935456*charge1*charge2/r
  </EnergyTerm>""" % numLJTypes
print ' <Function name="sigma" min="0" max="%d">' % (numLJTypes*numLJTypes-1)
for i in range(numLJTypes):
    for j in range(numLJTypes):
        print sigma[i][j] if j in sigma[i] else 0.5*(sigma[i][i]+sigma[j][j]),
    print
print ' </Function>'
for type in biotypes:
    if type[0] in chargeGroups:
        chargeGroup = chargeGroups[type[0]]
    else:
        chargeGroup = nextChargeGroupId
        nextChargeGroupId += 1
    fos = 0
    if type[0] in solvationParams:
        fos = solvationParams[type[0]]*4.184
    solvGroup = -1
    if type[0] in solvationGroups:
        solvGroup = solvationGroups[type[0]]
    if type[0] in volumeParams:
        values = volumeParams[type[0]]
        print '  <Atom type="%s" charge="%s" ljType="%s" epsilon="%g" radius="%g" chargeGroup="%d" solvGroup="%d" maxSav="%.3g" fos="%g" volReduction="%g"/>' % (type[0], charge[type[4]], int(type[3])-1, epsilon[type[3]], radius[type[3]], chargeGroup, solvGroup, values[1]/values[0], fos, values[2]/values[0])
    else:
        print '  <Atom type="%s" charge="%s" ljType="%s" epsilon="%g" radius="%g" chargeGroup="%d" solvGroup="-1" maxSav="1" fos="0" volReduction="1"/>' % (type[0], charge[type[4]], int(type[3])-1, epsilon[type[3]], radius[type[3]], chargeGroup)
print ' </CustomGBForce>'
print """ <Script>
import simtk.openmm as mm
gb = [f for f in [sys.getForce(i) for i in range(sys.getNumForces())] if type(f) == mm.CustomGBForce][0]

# Add Lennard-Jones exceptions.
for bond in data.bonds:
    gb.addExclusion(bond.atom1, bond.atom2)
for angle in data.angles:
    gb.addExclusion(angle[0], angle[2])

# Identify charge groups.
numChargeGroups = 0
chargeGroups = {}
groupParam = [gb.getParticleParameters(i)[4] for i in range(sys.getNumParticles())]
for atom in data.atoms:
    group = (atom.residue.index, groupParam[atom.index])
    if group in chargeGroups:
        groupId = chargeGroups[group]
    else:
        groupId = numChargeGroups
        numChargeGroups += 1
        chargeGroups[group] = groupId
    groupParam[atom.index] = groupId
for i in range(sys.getNumParticles()):
    params = list(gb.getParticleParameters(i))
    params[4] = groupParam[i]
    gb.setParticleParameters(i, params)

# Work out exclusions between charge groups.
groupExclusions = [set() for i in range(numChargeGroups)]
for angle in data.angles:
    for i in range(3):
        group1 = groupParam[angle[i]]
        for j in range(3):
            group2 = groupParam[angle[j]]
            groupExclusions[group1].add(group2)
includeFunction = [1.0 for i in range(numChargeGroups*numChargeGroups)]
for i in range(numChargeGroups):
    for j in range(numChargeGroups):
        if j in groupExclusions[i]:
            includeFunction[i*numChargeGroups+j] = 0.0
gb.addFunction('include', includeFunction, 0, numChargeGroups*numChargeGroups-1)
gb.setGlobalParameterDefaultValue(0, numChargeGroups)

# Assign backbone solvation groups.
numSolvGroups = 0
solvGroups = {}
groupParam = [gb.getParticleParameters(i)[5] for i in range(sys.getNumParticles())]
newParam = [-1]*len(groupParam)
residues = list(topology.residues())
isBackbone = [False for atom in data.atoms]
for i in range(len(residues)-1):
    groupAtoms = [atom for atom in residues[i].atoms() if groupParam[atom.index] == 0] + [atom for atom in residues[i+1].atoms() if groupParam[atom.index] == 1]
    groupId = numSolvGroups
    numSolvGroups += 1
    for atom in groupAtoms:
        newParam[atom.index] = groupId
        isBackbone[atom.index] = True

# Assign other solvation groups.
for atom in data.atoms:
    if isBackbone[atom.index] or groupParam[atom.index] == -1:
        continue
    group = (atom.residue.index, groupParam[atom.index])
    if group in solvGroups:
        groupId = solvGroups[group]
    else:
        groupId = numSolvGroups
        numSolvGroups += 1
        solvGroups[group] = groupId
    newParam[atom.index] = groupId
for i in range(sys.getNumParticles()):
    params = list(gb.getParticleParameters(i))
    params[5] = newParam[i]
    gb.setParticleParameters(i, params)
 </Script>"""
print '</ForceField>'
