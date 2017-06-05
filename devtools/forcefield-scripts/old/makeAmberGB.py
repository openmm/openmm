import sys
import xml.etree.ElementTree as etree

gbvalues = {'CT':(0.19,0.72),
            'CX':(0.19,0.72),
            'CI':(0.19,0.72),
            'C':(0.1875,0.72),
            'CA':(0.1875,0.72),
            'CM':(0.1875,0.72),
            'CS':(0.1875,0.72),
            'C4':(0.1875,0.72),
            'CC':(0.1875,0.72),
            'CV':(0.1875,0.72),
            'CW':(0.1875,0.72),
            'CR':(0.1875,0.72),
            'CB':(0.1875,0.72),
            'C*':(0.1875,0.72),
            'CN':(0.1875,0.72),
            'CK':(0.1875,0.72),
            'CP':(0.1875,0.72),
            'C5':(0.1875,0.72),
            'CQ':(0.1875,0.72),
            'N':(0.1706,0.79),
            'NA':(0.1706,0.79),
            'NB':(0.1706,0.79),
            'NC':(0.1706,0.79),
            'N*':(0.1706,0.79),
            'N2':(0.1706,0.79),
            'N3':(0.1625,0.79),
            'OW':(0.1535,0.85),
            'OH':(0.1535,0.85),
            'OS':(0.1535,0.85),
            'O':(0.148,0.85),
            'O2':(0.148,0.85),
            'S':(0.1775,0.96),
            'SH':(0.1775,0.96),
            'H':(0.115,0.85),
            'HW':(0.105,0.85),
            'HO':(0.105,0.85),
            'HS':(0.125,0.85),
            'HA':(0.125,0.85),
            'HC':(0.125,0.85),
            'H0':(0.125,0.85),
            'H1':(0.125,0.85),
            'H2':(0.125,0.85),
            'H3':(0.125,0.85),
            'HP':(0.125,0.85),
            'H4':(0.125,0.85),
            'H5':(0.125,0.85)}

tree = etree.parse(sys.argv[1])
typeMap = {}
for type in tree.getroot().find('AtomTypes').findall('Type'):
    typeMap[type.attrib['name']] = type.attrib['class']
print("<ForceField>")
print(" <GBSAOBCForce>")
for atom in tree.getroot().find('NonbondedForce').findall('Atom'):
    type = atom.attrib['type']
    if type in typeMap:
        atomClass = typeMap[type]
        if atomClass in gbvalues:
            values = gbvalues[atomClass]
            print("""   <Atom type="%s" charge="%s" radius="%g" scale="%g"/>""" % (type, atom.attrib['charge'], values[0], values[1]))
print(" </GBSAOBCForce>")
print("</ForceField>")

