#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
processTinkerForceField.py

   Convert TINKER force field files into xml files for use by pyopenmm

   (1) read residue template file
   (2) read TINKER parameter file
   (3) assign biotypes to each atom in residue template file
   (4) output force-field parameter file

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import xml.etree.ElementTree as etree
import sys
import shlex
import math
import datetime
import os.path

#=============================================================================================
# Ion list
#=============================================================================================

# biotype    2003    NA      "Sodium Ion"                      250
# biotype    2004    K       "Potassium Ion"                   251
# biotype    2005    MG      "Magnesium Ion"                   255
# biotype    2006    CA      "Calcium Ion"                     256
# biotype    2007    CL      "Chloride Ion"                    258

ions    = { 'Li+' :  ['LI', 249],
            'Na+' :  ['NA', 250],
            'K+'  :  ['K',  251],
            'Rb+' :  ['RB', 252],  
            'Cs+' :  ['CS', 253], 
            'Be+' :  ['BE', 254],
            'Mg+' :  ['MG', 255],
            'Ca+' :  ['CA', 256],
            'Zn+' :  ['ZN', 257],
            'Cl-' :  ['Cl', 258]
           }

atomTypes                    = {}
bioTypes                     = {}

#=============================================================================================
# Default 'constructor' for atoms 
#=============================================================================================

def getDefaultAtom( ):
 
    atom                     = dict();
    atom['tinkerLookupName'] = 'XXX'
    atom['type']             = -1
    atom['bonds']            = dict()

    return atom

#=============================================================================================
# Add bond to atomDict[]; atoms are added to atomDict[] if missing
#=============================================================================================

def addBond( atomDict, atom1, atom2 ):
    if( atom1 not in atomDict ):
        atomDict[atom1] = getDefaultAtom()
        
    if( atom2 not in atomDict ):
        atomDict[atom2] = getDefaultAtom()

    atomDict[atom2]['bonds'][atom1] = 1
    atomDict[atom1]['bonds'][atom2] = 1

#=============================================================================================
# Get atom dictionary from xml atom list
#=============================================================================================

def getXmlAtoms( atoms ):
 
    atomInfo = dict();
    for atom in atoms:
        name                               = atom.attrib['name']
        atomInfo[name]                     = getDefaultAtom()
        atomInfo[name]['tinkerLookupName'] = atom.attrib['tinkerLookupName']

    return atomInfo

#=============================================================================================
# Get bond dictionary from xml bond list
#=============================================================================================

def getXmlBonds( bonds ):
 
    bondInfo = dict();
    for bond in bonds:
        atom1                     = bond.attrib['from']
        atom2                     = bond.attrib['to']
        if( atom1 not in bondInfo ):
            bondInfo[atom1]       = dict()
        if( atom2 not in bondInfo ):
            bondInfo[atom2]       = dict()
    
        bondInfo[atom1][atom2] = 1
        bondInfo[atom2][atom1] = 1

    return bondInfo

#=============================================================================================
# Build entry for protein residue
#=============================================================================================

def buildProteinResidue( residueDict, atoms, bondInfo, abbr, loc, tinkerLookupName, include, residueName, type ):
    
    # residueDict[abbr]                         abbr=ALA, CALA, NALA, ...
    # residueDict[abbr]['atoms']                list if atom dict()
    # residueDict[abbr]['type']                 molecule type ('protein', 'nucleic acid', ...)
    # residueDict[abbr]['tinkerLookupName']     Tinker lookup name
    # residueDict[abbr]['residueName']          residueName
    # residueDict[abbr]['include']              include in output

    residueDict[abbr]                         = dict()
    residueDict[abbr]['atoms']                = atoms
    residueDict[abbr]['type']                 = type
    residueDict[abbr]['loc']                  = loc
    residueDict[abbr]['tinkerLookupName']     = tinkerLookupName
    residueDict[abbr]['residueName']          = residueName
    residueDict[abbr]['include']              = include

    # for each bond, add entry to 
    #   residueDict[abbr]['atoms'][atom]['bonds']
    #   residueDict[abbr]['atoms'][bondedAtom]['bonds']

    for atom in bondInfo:
        if( atom in residueDict[abbr]['atoms'] ):

            if( 'bonds' not in  residueDict[abbr]['atoms'][atom] ):
                residueDict[abbr]['atoms'][atom]['bonds'] = dict() 

            for bondedAtom in bondInfo[atom]:
                if( bondedAtom in  residueDict[abbr]['atoms'] ):
                    if( 'bonds' not in  residueDict[abbr]['atoms'][bondedAtom] ):
                        residueDict[abbr]['atoms'][bondedAtom]['bonds'] = dict() 
                    residueDict[abbr]['atoms'][bondedAtom]['bonds'][atom] = 1
                    residueDict[abbr]['atoms'][atom]['bonds'][bondedAtom] = 1
                else:
                    print "Error: bonded atom=%s not in residue=%s" % ( atom, abbr )
        else:
            print "Error: bonded atom=%s nt in residue=%s" % ( atom, abbr )

    return

#=============================================================================================
# Copy a bond (dict() copy)
#=============================================================================================

def copyBonds( bonds ):
    bondCopy = dict()
    for key in bonds.keys():
        bondCopy[key] = bonds[key]
    return bondCopy

#=============================================================================================
# Copy a atom (dict() copy, including the 'bonds' list)
#=============================================================================================

def copyAtom( atom ):
    atomCopy = dict()
    for key in atom.keys():
        if( key != 'bonds' ):
            atomCopy[key]      = atom[key]
        else:
            atomCopy['bonds']  = copyBonds( atom[key] )
    return atomCopy

#=============================================================================================
# Copy a residue, including atom list
#=============================================================================================

def copyProteinResidue( residue ):
    
    residueCopy                        = dict()
    residueCopy['atoms']               = dict()
    residueCopy['type']                = residue['type']
    residueCopy['loc']                 = residue['loc']
    residueCopy['tinkerLookupName']    = residue['tinkerLookupName']
    residueCopy['residueName']         = residue['residueName']
    residueCopy['include']             = residue['include']

    for atom in residue['atoms']:
        residueCopy['atoms'][atom] = copyAtom( residue['atoms'][atom] )

    return residueCopy

#=============================================================================================
# Build residue hash based on xml file
#=============================================================================================

def buildResidueDict( residueXmlFileName ):

    residueTree = etree.parse(residueXmlFileName)
    print "Read %s" % (residueXmlFileName)
    root        = residueTree.getroot()
    residueDict = dict()

    # residueDict[residueName] = dict()
    #      ['loc']
    #      ['type']
    #      ['atoms'] = dict()
    #          [atomName]    = dict()
    #              ['bonds'] = dict{}
    for residue in root.findall('Residue'):

        #  <Residue abbreviation="MET" loc="middle" type="protein" tinkerLookupName="Methionine" fullName="Methionine">
        abbr        = residue.attrib['abbreviation']
        loc         = residue.attrib['loc']
        type        = residue.attrib['type']
        tinkerName  = residue.attrib['tinkerLookupName']
        residueName = residue.attrib['fullName']
        isProtein   = False
        isWater     = False
        isDNA       = False
        isRNA       = False
        if type == 'protein':
            isProtein = True
        elif type == 'AmoebaWater':
            isWater = True
        elif type == 'dna':
            isDNA = True
        elif type == 'rna':
            isRNA = True
        else:
            continue

        atoms       = getXmlAtoms( residue.findall('Atom') ) 
        bondInfo    = getXmlBonds( residue.findall('Bond') ) 

        # if residue is an amino acid, then create CALA and NALA residues, in addition to non-termianal residue, and include approriate atoms
        # HXT is excluded from all residues

        if( isWater ):

            buildProteinResidue( residueDict, atoms, bondInfo, abbr, 'x', tinkerName, 1, 'HOH', 'water' )

        elif( isProtein ):

            buildProteinResidue( residueDict, atoms, bondInfo, abbr, 'm', tinkerName, 1, residueName, 'protein' )

            cResidueName                                                      = 'C' + abbr
            residueDict[cResidueName]                                         = copyProteinResidue( residueDict[abbr] )
            residueDict[cResidueName]['loc']                                  = 'c'
            if( residueDict[abbr]['tinkerLookupName'].find('(') > -1 ):
                begin = residueDict[abbr]['tinkerLookupName'].find('(')
                end   = residueDict[abbr]['tinkerLookupName'].find(')') + 1
                sub   = residueDict[abbr]['tinkerLookupName'][begin:end]
                if( sub == '(HD)' or sub == '(HE)' ):
                    residueDict[cResidueName]['tinkerLookupName']                 = 'C-Terminal ' + 'HIS ' + sub
                else:
                    residueDict[cResidueName]['tinkerLookupName']                 = 'C-Terminal ' + abbr + ' ' + sub
                print "tinkerLookupName %s %s" % ( abbr, residueDict[cResidueName]['tinkerLookupName'])
            else:
                residueDict[cResidueName]['tinkerLookupName']                 = 'C-Terminal ' + abbr
            residueDict[cResidueName]['atoms']['OXT']                         = copyAtom( residueDict[abbr]['atoms']['O'] )
            residueDict[cResidueName]['atoms']['OXT']['tinkerLookupName']     = 'OXT'
            residueDict[cResidueName]['atoms']['O']['tinkerLookupName']       = 'OXT'
            residueDict[cResidueName]['parent']                               = residueDict[abbr]
 
            nResidueName                                                      = 'N' + abbr
            residueDict[nResidueName]                                         = copyProteinResidue( residueDict[abbr] )
            residueDict[nResidueName]['loc']                                  = 'n'
            residueDict[nResidueName]['tinkerLookupName']                     = 'N-Terminal ' + abbr
            residueDict[nResidueName]['parent']                               = residueDict[abbr]

            if( abbr == 'PRO' ):
               #<Atom name="H" tinkerLookupName="HN" bonds="1" />

               residueDict[nResidueName]['atoms']['H2']                     = getDefaultAtom()
               residueDict[nResidueName]['atoms']['H3']                     = getDefaultAtom()

               residueDict[nResidueName]['atoms']['H2']['tinkerLookupName'] = 'HN'
               residueDict[nResidueName]['atoms']['H3']['tinkerLookupName'] = 'HN'

               addBond( residueDict[nResidueName]['atoms'], 'H2', 'N' )
               addBond( residueDict[nResidueName]['atoms'], 'H3', 'N' )
 
            else:
               residueDict[nResidueName]['atoms']['H2']     = copyAtom( residueDict[abbr]['atoms']['H'] )
               residueDict[nResidueName]['atoms']['H3']     = copyAtom( residueDict[abbr]['atoms']['H'] )

        elif isDNA or isRNA:

            buildProteinResidue( residueDict, atoms, bondInfo, abbr, loc, tinkerName, 1, residueName, type )

    print "Start Lookup XML FFFFinal\n\n"
    printXml = 1
    if( printXml ):
        print "<Residues>"
        for resName in sorted( residueDict.keys() ):
            if( 'include' in residueDict[resName] and residueDict[resName]['include'] ):
                type         = residueDict[resName]['type']
                loc          = residueDict[resName]['loc']
                tinkerLookupName   = residueDict[resName]['tinkerLookupName']
                fullName     = residueDict[resName]['residueName']
                outputString = """  <Residue abbreviation="%s" loc="%s" type="%s" tinkerLookupName="%s" fullName="%s">""" % (resName, loc, type, tinkerLookupName, fullName )
                print "%s" % outputString

                atomsInfo    = residueDict[resName]['atoms']
                for atomName in sorted( atomsInfo.keys() ):
                    tinkerLookupName = atomsInfo[atomName]['tinkerLookupName']
                    outputString = """  <Atom name="%s" tinkerLookupName="%s" />""" % (atomName, tinkerLookupName)
                    print "%s" % outputString

                includedBonds = dict()
                for atomName in sorted( atomsInfo.keys() ):
                    bondsInfo    = atomsInfo[atomName]['bonds']
                    for bondedAtom  in bondsInfo:
                        if( bondedAtom not in includedBonds or atomName not in includedBonds[bondedAtom] ):
                            outputString = """  <Bond from="%s" to="%s" />""" % (atomName, bondedAtom)
                            if( atomName not in includedBonds ):
                                includedBonds[atomName] = dict()
                            if( bondedAtom not in includedBonds ):
                                includedBonds[bondedAtom] = dict()
                            includedBonds[atomName][bondedAtom] = 1
                            includedBonds[bondedAtom][atomName] = 1
                            print "%s" % outputString
                print "</Residue>"
        print "</Residues>"

    return residueDict

#=============================================================================================
# Set biotype for each atom in residueDict
#=============================================================================================

def setBioTypes( bioTypes, residueDict ):

    for resname, res in residueDict.items():
        for atom in res['atoms']:
            atomLookup = res['atoms'][atom]['tinkerLookupName']
            resLookup = []
            if res['type'] == 'dna':
                if res['loc'] in ('5', 'N'):
                    resLookup.append("5'-Hydroxyl DNA")
                if res['loc'] in ('3', 'N'):
                    resLookup.append("3'-Hydroxyl DNA")
                resLookup.append("Phosphodiester DNA")
            if res['type'] == 'rna':
                if res['loc'] in ('5', 'N'):
                    resLookup.append("5'-Hydroxyl RNA")
                if res['loc'] in ('3', 'N'):
                    resLookup.append("3'-Hydroxyl RNA")
                resLookup.append("Phosphodiester RNA")
            resLookup.append(res['tinkerLookupName'])
            for suffix in resLookup:
                lookupName = atomLookup+'_'+suffix
                if lookupName in bioTypes:
                    break
            if lookupName in bioTypes:
                res['atoms'][atom]['type'] = bioTypes[lookupName][3]
            else:
                print "For %s lookupName=%s not in biotype" % (atom,lookupName)
                if( 'parent' in res ):
                    lookupName =  res['atoms'][atom]['tinkerLookupName']  + '_' +  res['parent']['tinkerLookupName']
                    if( lookupName in bioTypes ):
                        res['atoms'][atom]['type'] = bioTypes[lookupName][3]
                    else:
                        print "Missing lookupName=%s from biotype" % (lookupName)
    return 0

#=============================================================================================
# Add multipole for forces[]; added entry is a list of axis info [kz, kx, ky] and another 
# list of multipoles [charge, dipole, quadrupole]
#=============================================================================================

def addMultipole( lineIndex, allLines, forces ):

    if( 'multipole' not in forces ):
        forces['multipole'] = []

    # axis indices and charge

    fields                  = allLines[lineIndex]
    multipoles              = [ fields[-1] ]
    axisInfo                = fields[1:-1]

    # dipole

    lineIndex              += 1
    fields                  = allLines[lineIndex]
    multipoles.append( fields[0] )
    multipoles.append( fields[1] )
    multipoles.append( fields[2] )

    # quadrupole

    lineIndex              += 1
    fields                  = allLines[lineIndex]
    multipoles.append( fields[0] )

    lineIndex              += 1
    fields                  = allLines[lineIndex]
    multipoles.append( fields[0] )
    multipoles.append( fields[1] )

    lineIndex              += 1
    fields                  = allLines[lineIndex]
    multipoles.append( fields[0] )
    multipoles.append( fields[1] )
    multipoles.append( fields[2] )

    lineIndex              += 1

    # save info

    multipoleInfo           = [ axisInfo, multipoles ]
    forces['multipole'].append( multipoleInfo )

    return (lineIndex)

#=============================================================================================
# Add tortor parameters/grid to forces[]; format of each entry is [ first tortor line, grid ]
#=============================================================================================

def addTorTor( lineIndex, allLines, forces ):

    if( 'tortors' not in forces ):
        forces['tortors'] = []

    fields         = allLines[lineIndex]
    tortorInfo     = fields[1:]

    # read grid lines

    lastGridLine   = lineIndex + int(fields[6])*int(fields[7])
    grid           = []
    while( lineIndex < lastGridLine ):
        lineIndex += 1
        grid.append( allLines[lineIndex] )

    forces['tortors'].append( [ tortorInfo, grid ] )

    return (lineIndex)

#=============================================================================================

residueXmlFileName = 'residuesFinal.xml'
residueDict        = buildResidueDict( residueXmlFileName )
    
#=============================================================================================

# recognizedForces[] contain raw list entries from TINKER parameter file

resAtomTypes                             = {}
forces                                   = {}
recognizedForces                         = {}
recognizedForces['bond']                 = 1
recognizedForces['angle']                = 1
recognizedForces['strbnd']               = 1
recognizedForces['ureybrad']             = 1
recognizedForces['opbend']               = 1
recognizedForces['torsion']              = 1
recognizedForces['pitors']               = 1
recognizedForces['vdw']                  = 1
recognizedForces['polarize']             = 1
recognizedForces['tortors']              = addTorTor
recognizedForces['multipole']            = addMultipole

#=============================================================================================

# recognizedScalars[] contain raw scalar entries from TINKER parameter file

scalars                                  = {}
recognizedScalars                        = {}
recognizedScalars['forcefield']          = '-2.55'
recognizedScalars['bond-cubic']          = '-2.55'
recognizedScalars['bond-quartic']        = '3.793125'
recognizedScalars['angle-cubic']         = '-0.014'
recognizedScalars['angle-quartic']       = '0.000056'
recognizedScalars['angle-pentic']        = '-0.0000007'
recognizedScalars['angle-sextic']        = '0.000000022'
recognizedScalars['opbendtype']          = 'ALLINGER'
recognizedScalars['opbend-cubic']        = '-0.014'
recognizedScalars['opbend-quartic']      = '0.000056'
recognizedScalars['opbend-pentic']       = '-0.0000007'
recognizedScalars['opbend-sextic']       = '0.000000022'
recognizedScalars['torsionunit']         = '0.5'
recognizedScalars['vdwtype']             = 'BUFFERED-14-7'
recognizedScalars['radiusrule']          = 'CUBIC-MEAN'
recognizedScalars['radiustype']          = 'R-MIN'
recognizedScalars['radiussize']          = 'DIAMETER'
recognizedScalars['epsilonrule']         = 'HHG'
recognizedScalars['dielectric']          = '1.0'
recognizedScalars['polarization']        = 'MUTUAL'
recognizedScalars['vdw-13-scale']        = '0.0'
recognizedScalars['vdw-14-scale']        = '1.0'
recognizedScalars['vdw-15-scale']        = '1.0'
recognizedScalars['mpole-12-scale']      = '0.0'
recognizedScalars['mpole-13-scale']      = '0.0'
recognizedScalars['mpole-14-scale']      = '0.4'
recognizedScalars['mpole-15-scale']      = '0.8'
recognizedScalars['polar-12-scale']      = '0.0'
recognizedScalars['polar-13-scale']      = '0.0'
recognizedScalars['polar-14-scale']      = '1.0'
recognizedScalars['polar-15-scale']      = '1.0'
recognizedScalars['polar-14-intra']      = '0.5'
recognizedScalars['direct-11-scale']     = '0.0'
recognizedScalars['direct-12-scale']     = '1.0'
recognizedScalars['direct-13-scale']     = '1.0'
recognizedScalars['direct-14-scale']     = '1.0'
recognizedScalars['mutual-11-scale']     = '1.0'
recognizedScalars['mutual-12-scale']     = '1.0'
recognizedScalars['mutual-13-scale']     = '1.0'
recognizedScalars['mutual-14-scale']     = '1.0'

#=============================================================================================
# get all 'interesting' lines in file

allLines                                 = []
for line in open(sys.argv[1]):
    try:
        fields = shlex.split(line)
    except:
        continue
    if len(fields) == 0:
        continue
    if fields[0][0] == '#':
        continue
    allLines.append( fields )

#=============================================================================================

# load lines in lists/scalar values

lineIndex = 0
while lineIndex < len( allLines ):

    fields = allLines[lineIndex]

    if fields[0] == 'atom':
        if( fields[3] in ions ):
            ionInfo = ions[fields[3]]
            element = ionInfo[0]
            ionInfo[1] = int(fields[1])
        else:
            element = fields[3][0]
        atomTypes[int(fields[1])] = (fields[2], element, fields[6])
        lineIndex += 1

    elif fields[0] == 'biotype':

        lookUp            = fields[2] + '_' + fields[3]
        if lookUp in bioTypes:
            # Workaround for Tinker using the same name but different types for H2', H2'', and for H5', H5''
            lookUp = fields[2]+'*_'+fields[3]
        bioTypes[lookUp]  = fields[1:]
        lineIndex        += 1

    elif fields[0] in recognizedForces:
        if( recognizedForces[fields[0]] == 1 ):
            if( fields[0] not in forces ):
                forces[fields[0]] = []
            forces[fields[0]].append( fields[1:] )
            lineIndex += 1
        else:
            lineIndex = recognizedForces[fields[0]]( lineIndex, allLines, forces )

    elif fields[0] in recognizedScalars:
        scalars[fields[0]] = fields[1]
        lineIndex += 1
    else:
        print "Field %s not recognized: line=<%s>" % ( fields[0], allLines[lineIndex] )
        lineIndex += 1

#=============================================================================================

# set biotypes for all atoms

setBioTypes( bioTypes, residueDict )

#=============================================================================================

# open force field xml file for output

tinkerXmlFileName           = scalars['forcefield']
tinkerXmlFileName          += '.xml'
tinkerXmlFile               = open( tinkerXmlFileName, 'w' )
print "Opened %s." % (tinkerXmlFileName)

gkXmlFileName              = scalars['forcefield']
gkXmlFileName             += '_gk.xml'
gkXmlFile                  = open( gkXmlFileName, 'w' )
print "Opened %s." % (gkXmlFileName)

today = datetime.date.today().isoformat()
sourceFile = os.path.basename(sys.argv[1])
header = """ <Info>
  <Source>%s</Source>
  <DateGenerated>%s</DateGenerated>
  <Reference></Reference>
 </Info>
""" % (sourceFile, today)

gkXmlFile.write( "<ForceField>\n" )
gkXmlFile.write(header)
tinkerXmlFile.write( "<ForceField>\n" )
tinkerXmlFile.write(header)
tinkerXmlFile.write( " <AtomTypes>\n")

if( scalars['forcefield'].find( 'AMOEBA' ) > -1 ):
    isAmoeba = 1
else:
    isAmoeba = 0

#=============================================================================================

# atom type/class

#         atmType  class   name  name                    atomicNo.     mass valence     
# atom          1    1    N     "Glycine N"                    7    14.003    3
# atom          2    2    CA    "Glycine CA"                   6    12.000    4
# atom          3    3    C     "Glycine C"                    6    12.000    3
# atom          4    4    HN    "Glycine HN"                   1     1.008    1
# atom          5    5    O     "Glycine O"                    8    15.995    1

# atom        380   73    O     "AMOEBA Water O"               8    15.999    2
# atom        381   74    H     "AMOEBA Water H"               1     1.008    1
# atom        383   76    Na+   "Sodium Ion Na+"              11    22.990    0
# atom        384   77    K+    "Potassium Ion K+"            19    39.098    0
# atom        385   78    Rb+   "Rubidium Ion Rb+"            37    85.468    0
# atom        386   79    Cs+   "Cesium Ion Cs+"              55   132.905    0
# atom        387   80    Be+   "Beryllium Ion Be+2"           4     9.012    0
# atom        388   81    Mg+   "Magnesium Ion Mg+2"          12    24.305    0
# atom        389   82    Ca+   "Calcium Ion Ca+2"            20    40.078    0
# atom        390   83    Cl-   "Chloride Ion Cl-"            17    35.453    0
  

#            biotype                                          atmType
# biotype       1    N       "Glycine"                           1     
# biotype       2    CA      "Glycine"                           2     
# biotype       3    C       "Glycine"                           3     
# biotype       4    HN      "Glycine"                           4     
# biotype       5    O       "Glycine"                           5     

# biotype    2001    O       "Water"                           380
# biotype    2002    H       "Water"                           381
# biotype    2003    NA      "Sodium Ion"                      383   
# biotype    2004    K       "Potassium Ion"                   384   
# biotype    2005    MG      "Magnesium Ion"                   388   
# biotype    2006    CA      "Calcium Ion"                     389   
# biotype    2007    CL      "Chloride Ion"                    390   
  
if( isAmoeba ):
    for type in sorted(atomTypes):
        outputString = """  <Type name="%s" class="%s" element="%s" mass="%s"/>""" % (type, atomTypes[type][0], atomTypes[type][1], atomTypes[type][2])
        tinkerXmlFile.write( "%s\n" % (outputString) )
else:
    for type in sorted(atomTypes):
        outputString = """  <Type name="%s" class="%s" mass="%s"/>""" % (type, atomTypes[type][0], atomTypes[type][1])
        tinkerXmlFile.write( "%s\n" % (outputString) )

tinkerXmlFile.write( " </AtomTypes>\n")

#=============================================================================================

# residues

tinkerXmlFile.write( " <Residues>\n" )
for resname, res in sorted(residueDict.items()):
    if res['include']:
        outputString = """  <Residue name="%s">""" % (resname)
        tinkerXmlFile.write( "%s\n" % (outputString) )
        atomIndex    = dict()
        atomCount    = 0
        for atom in sorted( res['atoms'].keys() ):
            type  = res['atoms'][atom]['type']
            typeI = int( type )
            if( typeI < 0 ):
                print "Error: type=%s for atom=%s of residue=%s" % (type, atom, resname)
            tag  = "   <Atom name=\"%s\" type=\"%s\" />" % (atom, type)
            atomIndex[atom]  = atomCount
            atomCount       += 1
            tinkerXmlFile.write( "%s\n" % (tag) )
    
        includedBonds = dict()
        for atomName in sorted( res['atoms'].keys() ):
            bondsInfo    = res['atoms'][atomName]['bonds']
            for bondedAtom  in bondsInfo:
                if( bondedAtom not in includedBonds or atomName not in includedBonds[bondedAtom] ):
                    outputString = """   <Bond from="%s" to="%s" />""" % (str(atomIndex[atomName]), str(atomIndex[bondedAtom]))
                    if( atomName not in includedBonds ):
                        includedBonds[atomName] = dict()
                    if( bondedAtom not in includedBonds ):
                        includedBonds[bondedAtom] = dict()
                    includedBonds[atomName][bondedAtom] = 1
                    includedBonds[bondedAtom][atomName] = 1
                    tinkerXmlFile.write( "%s\n" % (outputString) )

        outputStrings = []
        if res['type'] in ('rna', 'dna'):
            if res['loc'] in ('middle', '3'):
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['P']) ))
            if res['loc'] in ('middle', '5'):
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex["O3'"]) ))
        else:
            if res['loc'] == 'm':
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['N']) ))
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['C']) ) )
            if res['loc'] == 'n':
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['C']) ) )
            if res['loc'] == 'c':
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['N']) ) )
            if resname.find('CYX') > -1:
                outputStrings.append( """   <ExternalBond from="%s" />""" % (str(atomIndex['SG']) ) )
        for outputString in outputStrings:
            tinkerXmlFile.write( "%s\n" % (outputString) )

        tinkerXmlFile.write( " </Residue>\n" )

# End caps

tinkerXmlFile.write('  <Residue name="ACE">\n')
tinkerXmlFile.write('   <Atom name="HH31" type="%s"/>\n' % bioTypes['H_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Atom name="CH3" type="%s"/>\n' % bioTypes['CH3_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Atom name="HH32" type="%s"/>\n' % bioTypes['H_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Atom name="HH33" type="%s"/>\n' % bioTypes['H_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Atom name="C" type="%s"/>\n' % bioTypes['C_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Atom name="O" type="%s"/>\n' % bioTypes['O_Acetyl N-Terminus'][3])
tinkerXmlFile.write('   <Bond from="0" to="1"/>\n')
tinkerXmlFile.write('   <Bond from="1" to="2"/>\n')
tinkerXmlFile.write('   <Bond from="1" to="3"/>\n')
tinkerXmlFile.write('   <Bond from="1" to="4"/>\n')
tinkerXmlFile.write('   <Bond from="4" to="5"/>\n')
tinkerXmlFile.write('   <ExternalBond from="4"/>\n')
tinkerXmlFile.write('  </Residue>\n')
tinkerXmlFile.write('  <Residue name="NME">\n')
tinkerXmlFile.write('   <Atom name="N" type="%s"/>\n' % bioTypes['N_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Atom name="H" type="%s"/>\n' % bioTypes['HN_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Atom name="CH3" type="%s"/>\n' % bioTypes['CH3_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Atom name="HH31" type="%s"/>\n' % bioTypes['H_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Atom name="HH32" type="%s"/>\n' % bioTypes['H_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Atom name="HH33" type="%s"/>\n' % bioTypes['H_N-MeAmide C-Terminus'][3])
tinkerXmlFile.write('   <Bond from="0" to="1"/>\n')
tinkerXmlFile.write('   <Bond from="0" to="2"/>\n')
tinkerXmlFile.write('   <Bond from="2" to="3"/>\n')
tinkerXmlFile.write('   <Bond from="2" to="4"/>\n')
tinkerXmlFile.write('   <Bond from="2" to="5"/>\n')
tinkerXmlFile.write('   <ExternalBond from="0"/>\n')
tinkerXmlFile.write('  </Residue>\n')


# ions
 
for ion,ionInfo in ions.iteritems():
    outputString  = """  <Residue name="%s">\n""" % (ionInfo[0])
    outputString += """   <Atom name="%s" type="%s"/>\n""" % (ionInfo[0], str(ionInfo[1]))
    outputString += """  </Residue>\n""" 
    tinkerXmlFile.write( "%s" % (outputString) )
    
tinkerXmlFile.write( " </Residues>\n" )

radian              = 57.2957795130
if( isAmoeba ):

#=============================================================================================

    # AmoebaBondForce

    cubic           = 10.*float(scalars['bond-cubic']) 
    quartic         = 100.*float(scalars['bond-quartic']) 
    outputString    = """ <AmoebaBondForce bond-cubic="%s" bond-quartic="%s">""" %( str(cubic), str(quartic) )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    bonds        = forces['bond']
    for bond in bonds:
       length       = float(bond[3])*0.1
       k            = float(bond[2])*100.0*4.184
       outputString = """  <Bond class1="%s" class2="%s" length="%s" k="%s"/>""" % (bond[0], bond[1], str(length), str(k))
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaBondForce>\n" )

#=============================================================================================

    # AmoebaAngleForce

    cubic           = float(scalars['angle-cubic']) 
    quartic         = float(scalars['angle-quartic']) 
    pentic          = float(scalars['angle-pentic']) 
    sextic          = float(scalars['angle-sextic']) 
    outputString    = """ <AmoebaAngleForce angle-cubic="%s" angle-quartic="%s" angle-pentic="%s" angle-sextic="%s">""" %( str(cubic), str(quartic), str(pentic), str(sextic) )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    angles          = forces['angle']
    radian          = 57.2957795130
    radian2         = 4.184/(radian*radian)
    for angle in angles:
       k            = float(angle[3])*radian2
       outputString = """  <Angle class1="%s" class2="%s" class3="%s" k="%s" angle1="%s" """ % (angle[0], angle[1], angle[2], str(k), angle[4] ) 
       if( len(angle) > 5 ):
           outputString += """  angle2="%s" """ % (angle[5])

       if( len(angle) > 6 ):
           outputString += """  angle3="%s" """ % (angle[6])
       outputString += " /> "

       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaAngleForce>\n" )

#=============================================================================================

    # AmoebaOutOfPlaneBendForce

    cubic           = float(scalars['opbend-cubic']) 
    quartic         = float(scalars['opbend-quartic']) 
    pentic          = float(scalars['opbend-pentic']) 
    sextic          = float(scalars['opbend-sextic']) 
    type            = scalars['opbendtype'] 
    outputString    = """ <AmoebaOutOfPlaneBendForce type="%s" opbend-cubic="%s" opbend-quartic="%s" opbend-pentic="%s" opbend-sextic="%s">""" % ( 
                            type, str(cubic), str(quartic), str(pentic), str(sextic) )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    opbends         = forces['opbend']
    radian2         = 4.184/(radian*radian)
    for opbend in opbends:
       k            = float(opbend[4])*radian2
       outputString = """  <Angle class1="%s" class2="%s" class3="%s" class4="%s" k="%s"/>""" % (opbend[0], opbend[1], opbend[2],  opbend[3], str(k))
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaOutOfPlaneBendForce>\n" )

#=============================================================================================

    # AmoebaTorsionForce

    torsionUnit      = float(scalars['torsionunit']) 
    outputString     = """ <PeriodicTorsionForce>"""
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    torsions         = forces['torsion']
    conversion       = 4.184*torsionUnit
    for torsion in torsions:
       outputString  = """  <Proper class1="%s" class2="%s" class3="%s" class4="%s" """ % (torsion[0], torsion[1], torsion[2],  torsion[3])
       startIndex    = 4
       for ii in range(0,3):
          torsionSuffix            = str(ii+1)
          amplitudeAttributeName   = 'k'+torsionSuffix
          angleAttributeName       = 'phase'+torsionSuffix
          periodicityAttributeName = 'periodicity'+torsionSuffix
          amplitude                = float(torsion[startIndex])*conversion
          angle                    = float(torsion[startIndex+1])/radian
          periodicity              = int(torsion[startIndex+2])
          outputString            += """  %s="%s" %s="%s" %s="%s" """ % (amplitudeAttributeName, str(amplitude), angleAttributeName, str(angle), periodicityAttributeName, str(periodicity))
          startIndex              += 3
       outputString += "/>"
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </PeriodicTorsionForce>\n" )

#=============================================================================================

    # AmoebaPiTorsionForce

    piTorsionUnit      = 1.0
    outputString       = """ <AmoebaPiTorsionForce piTorsionUnit="%s">""" % ( piTorsionUnit )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    piTorsions         = forces['pitors']
    conversion         = 4.184*piTorsionUnit
    for piTorsion in piTorsions:
       k             = float(piTorsion[2])*conversion
       outputString  = """  <PiTorsion class1="%s" class2="%s" k="%s" />""" % (piTorsion[0], piTorsion[1], str(k) )
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaPiTorsionForce>\n" )

#=============================================================================================

    # AmoebaStretchBendForce

    stretchBendUnit      = 1.0
    outputString         = """ <AmoebaStretchBendForce stretchBendUnit="%s">""" % ( stretchBendUnit )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    conversion           = 41.84/radian
    stretchBends         = forces['strbnd']
    for stretchBend in stretchBends:
       k1            = float(stretchBend[3])*conversion
       k2            = float(stretchBend[4])*conversion
       outputString  = """  <StretchBend class1="%s" class2="%s" class3="%s" k1="%s" k2="%s" />""" % (stretchBend[0], stretchBend[1], stretchBend[2], str(k1), str(k2) )
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( "</AmoebaStretchBendForce>\n" )

#=============================================================================================

    # AmoebaTorsionTorsionForce

    torsionTorsionUnit      = 1.0
    outputString            = """ <AmoebaTorsionTorsionForce >"""
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    torsionTorsions         = forces['tortors']
    for (index, torsionTorsion) in enumerate(torsionTorsions):
       torInfo       = torsionTorsion[0]
       grid          = torsionTorsion[1]
       outputString  = """  <TorsionTorsion class1="%s" class2="%s" class3="%s" class4="%s" class5="%s" grid="%s" nx="%s" ny="%s" />""" % (torInfo[0], torInfo[1], torInfo[2], torInfo[3], torInfo[4], str(index),
                              torInfo[5], torInfo[6] )
       tinkerXmlFile.write( "%s\n" % (outputString ) )

    for (index, torsionTorsion) in enumerate(torsionTorsions):
       torInfo       = torsionTorsion[0]
       grid          = torsionTorsion[1]
       outputString  = """  <TorsionTorsionGrid grid="%s" nx="%s" ny="%s" >""" % (str(index), torInfo[5], torInfo[6] )
       tinkerXmlFile.write( "%s\n" % (outputString ) )
       for (gridIndex, gridEntry) in enumerate(grid):
           print "Gxx %d  %s" % ( gridIndex, str(gridEntry) )
           if( len( gridEntry ) > 5 ):
               f   = float( gridEntry[2] )*4.184
               fx  = float( gridEntry[3] )*4.184
               fy  = float( gridEntry[4] )*4.184
               fxy = float( gridEntry[5] )*4.184
               outputString  = """  <Grid angle1="%s" angle2="%s" f="%s" fx="%s" fy="%s" fxy="%s" />""" % ( gridEntry[0], gridEntry[1], str(f), str(fx), str(fy), str(fxy) )
               tinkerXmlFile.write( "  %s\n" % (outputString ) )
           elif( len( gridEntry ) > 2 ):
               f   = float( gridEntry[2] )*4.184
               outputString  = """  <Grid angle1="%s" angle2="%s" f="%s" />""" % ( gridEntry[0], gridEntry[1], str(f) )
               tinkerXmlFile.write( "  %s\n" % (outputString ) )
       outputString  = '</TorsionTorsionGrid >'
       tinkerXmlFile.write( "%s\n" % (outputString ) )

    tinkerXmlFile.write( "</AmoebaTorsionTorsionForce>\n" )

#=============================================================================================

    # AmoebaVdwForce

    outputString         = """ <AmoebaVdwForce type="%s" radiusrule="%s" radiustype="%s" radiussize="%s" epsilonrule="%s" vdw-13-scale="%s" vdw-14-scale="%s" vdw-15-scale="%s" >""" % (
         scalars['vdwtype'], scalars['radiusrule'], scalars['radiustype'], scalars['radiussize'], scalars['epsilonrule'], scalars['vdw-13-scale'], scalars['vdw-14-scale'], scalars['vdw-15-scale'] )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    vdws                 = forces['vdw']
    for vdw in vdws:
       sigma             = float(vdw[1])*0.1
       epsilon           = float(vdw[2])*4.184
       if( len(vdw) > 3 ): 
           reduction = vdw[3]
       else:
           reduction = 1.0
       outputString      = """  <Vdw class="%s" sigma="%s" epsilon="%s" reduction="%s" /> """ % (vdw[0], str(sigma), str(epsilon), str(reduction) )
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaVdwForce>\n" )

#=============================================================================================

    # AmoebaMultipoleForce

    scalarList                  = dict()
    scalarList['mpole12Scale']  = recognizedScalars['mpole-12-scale']
    scalarList['mpole13Scale']  = recognizedScalars['mpole-13-scale']
    scalarList['mpole14Scale']  = recognizedScalars['mpole-14-scale']
    scalarList['mpole15Scale']  = recognizedScalars['mpole-15-scale']

    scalarList['polar12Scale']  = recognizedScalars['polar-12-scale']
    scalarList['polar13Scale']  = recognizedScalars['polar-13-scale']
    scalarList['polar14Scale']  = recognizedScalars['polar-14-scale']
    scalarList['polar15Scale']  = recognizedScalars['polar-15-scale']
    scalarList['polar14Intra']  = recognizedScalars['polar-14-intra']

    scalarList['direct11Scale'] = recognizedScalars['direct-11-scale']
    scalarList['direct12Scale'] = recognizedScalars['direct-12-scale']
    scalarList['direct13Scale'] = recognizedScalars['direct-13-scale']
    scalarList['direct14Scale'] = recognizedScalars['direct-14-scale']

    scalarList['mutual11Scale'] = recognizedScalars['mutual-11-scale']
    scalarList['mutual12Scale'] = recognizedScalars['mutual-12-scale']
    scalarList['mutual13Scale'] = recognizedScalars['mutual-13-scale']
    scalarList['mutual14Scale'] = recognizedScalars['mutual-14-scale']

    outputString         = """ <AmoebaMultipoleForce """
    for key in sorted( scalarList.keys() ): 
        outputString    += """ %s="%s" """ % ( key, scalarList[key] )
    outputString        += """ > """
    tinkerXmlFile.write( "%s\n" % (outputString ) )

    multipoleArray       = forces['multipole']
    bohr                 = 0.52917720859
    dipoleConversion     = 0.1*bohr
    quadrupoleConversion = 0.01*bohr*bohr/3.0
    for multipoleInfo in multipoleArray:
       axisInfo          = multipoleInfo[0]
       multipoles        =  multipoleInfo[1]
       outputString      = """  <Multipole type="%s" """ % (axisInfo[0] )
       axisInfoLen       = len(axisInfo)

       if( axisInfoLen > 1 ):
           outputString += """kz="%s" """ % ( axisInfo[1] )

       if( axisInfoLen > 2 ):
           outputString += """kx="%s" """ % ( axisInfo[2] )

       if( axisInfoLen > 3 ):
           outputString += """ky="%s" """ % ( axisInfo[3] )

       outputString += """c0="%s" d1="%s" d2="%s" d3="%s" q11="%s" q21="%s" q22="%s" q31="%s" q32="%s" q33="%s"  """ % ( multipoles[0],
                                    str(     dipoleConversion*float(multipoles[1]) ), str(    dipoleConversion*float(multipoles[2])), str(    dipoleConversion*float(multipoles[3])),
                                    str( quadrupoleConversion*float(multipoles[4]) ), str(quadrupoleConversion*float(multipoles[5])), str(quadrupoleConversion*float(multipoles[6])),
                                    str( quadrupoleConversion*float(multipoles[7]) ), str(quadrupoleConversion*float(multipoles[8])), str(quadrupoleConversion*float(multipoles[9])) )
       outputString     += "/>"
       tinkerXmlFile.write( "%s\n" % (outputString ) )

    polarizeArray = forces['polarize']
    
    polarityConversion = 0.001
    m = {}
    for polarize in polarizeArray:
       m[polarize[0]] = []
       outputString      = """  <Polarize type="%s" polarizability="%s" thole="%s" """ % (polarize[0], str(polarityConversion*float(polarize[1])), polarize[2] )
       for ii in range( 3, len(polarize) ):
          outputString  += """pgrp%d="%s" """ % (ii-2,polarize[ii])
          m[polarize[0]].append(polarize[ii])
          
       outputString     += "/>"
       tinkerXmlFile.write( "%s\n" % (outputString ) )
       print m[polarize[0]]
    for t in sorted(m):
        for k in m[t]:
            if t not in m[k]:
                print t, k

    tinkerXmlFile.write( " </AmoebaMultipoleForce>\n" )

#=============================================================================================

    # AmoebaGeneralizedKirkwoodForce

    solventDielectric = 78.3
    soluteDielectric  = 1.0
    includeCavityTerm = 1
    probeRadius       = 0.14
    surfaceAreaFactor = -6.0*3.1415926535*0.0216*1000.0*0.4184
    outputString      = """ <AmoebaGeneralizedKirkwoodForce solventDielectric="%s" soluteDielectric="%s" includeCavityTerm="%s" probeRadius="%s" surfaceAreaFactor="%s">""" % (
                           str(solventDielectric), str(soluteDielectric), str(includeCavityTerm), str(probeRadius), str(surfaceAreaFactor) )
    gkXmlFile.write( "%s\n" % (outputString ) )

    # radii are set in forcefield.py

    for type in sorted( atomTypes ):
        print "atom type=%s  %s" % ( str(type), str(atomTypes[type]) )

    for type in sorted( bioTypes ):
        print "bio type=%s  %s" % ( str(type), str(bioTypes[type]) )

    multipoleArray       = forces['multipole']
    for multipoleInfo in multipoleArray:
       axisInfo          = multipoleInfo[0]
       multipoles        =  multipoleInfo[1]
       type              = int(axisInfo[0])
       shct              = 0.8
       if( type in atomTypes ):
           element =  atomTypes[type][1]
           if( element == 'H' ):
               shct = 0.85
           elif( element == 'C' ):
               shct = 0.72
           elif( element == 'N' ):
               shct = 0.79
           elif( element == 'O' ):
               shct = 0.85
           elif( element == 'F' ):
               shct = 0.88
           elif( element == 'P' ):
               shct = 0.86
           elif( element == 'S' ):
               shct = 0.96
           elif( element == 'Fe' ):
               shct = 0.88
           else:
               print "Warning no overlap scale factor for type=%d element=%s" % (type, element)
       else:
           print "Warning no overlap scale factor for type=%d " % (type)

       outputString      = """  <GeneralizedKirkwood type="%s" charge="%s" shct="%s"  /> """ % ( axisInfo[0], multipoles[0],  str(shct) )
       gkXmlFile.write( "%s\n" % (outputString ) )
    gkXmlFile.write( " </AmoebaGeneralizedKirkwoodForce>\n" )

#=============================================================================================

    # AmoebaWcaDispersionForce

    epso         = 0.1100
    epsh         = 0.0135
    rmino        = 1.7025
    rminh        = 1.3275
    awater       = 0.033428
    slevy        = 1.0 
    dispoff      = 0.26
    shctd        = 0.81

    outputString         = """ <AmoebaWcaDispersionForce epso="%s" epsh="%s" rmino="%s" rminh="%s" awater="%s" slevy="%s"  dispoff="%s" shctd="%s" >""" % (
                           str(epso*4.184), str( epsh*4.184), str( rmino*0.1), str( rminh*0.1), str( 1000.0*awater), str( slevy), str( 0.1*dispoff), str( shctd ) )
    gkXmlFile.write( "%s\n" % (outputString ) )
    vdws                 = forces['vdw']
    convert              = 0.1
    if( scalars['radiustype'] == 'SIGMA' ):
        convert         *= 1.122462048309372

    if( scalars['radiussize'] == 'DIAMETER' ):
        convert         *= 0.5

    for vdw in vdws:
       sigma             = float(vdw[1])
       sigma            *= convert
       epsilon           = float(vdw[2])*4.184
       outputString      = """  <WcaDispersion class="%s" radius="%s" epsilon="%s" /> """ % ( vdw[0], str(sigma), str(epsilon) )
       gkXmlFile.write( "%s\n" % (outputString ) )
    gkXmlFile.write( " </AmoebaWcaDispersionForce>\n" )

#=============================================================================================

    # AmoebaUreyBradleyForce

    cubic        = 0.0
    quartic      = 0.0

    outputString         = """ <AmoebaUreyBradleyForce cubic="%s" quartic="%s"  >""" % ( str(cubic), str(quartic) )
    tinkerXmlFile.write( "%s\n" % (outputString ) )
    ubs                  = forces['ureybrad']
    for ub in ubs:
       k                 = float(ub[3])*4.184*100.0
       d                 = float(ub[4])*0.1
       outputString      = """  <UreyBradley class1="%s" class2="%s" class3="%s" k="%s" d="%s" /> """ % ( ub[0],  ub[1],  ub[2], str(k), str(d) )
       tinkerXmlFile.write( "%s\n" % (outputString ) )
    tinkerXmlFile.write( " </AmoebaUreyBradleyForce>\n" )

#=============================================================================================

tinkerXmlFile.write("</ForceField>\n")
gkXmlFile.write("</ForceField>\n")
tinkerXmlFile.close()
gkXmlFile.close()

