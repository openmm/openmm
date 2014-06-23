import numpy as np
import math, sys, string, shlex
import xml.etree.ElementTree as etree
from simtk.openmm.app import PDBFile

energy_conversion =  4.184 # convert kcal to kJ
length_conversion = 0.1 #convert angstrom to nm
angle_conversion = np.pi/180. #convert degree to rad

Root = "/opt/schrodinger/desmond-v31023/data/viparr/ff3/charmm22star/"

tinker_root = "/opt/pymd/KyleCode/App/Charmm22Star/TinkerFiles/"
tinker_impropers = np.loadtxt(tinker_root+"TinkerImpropers.dat")
tinker_atoms = np.loadtxt(tinker_root+"TinkerAtoms.dat","str")

def ProcessTinkerImpropers():
	atom_dict = {}
	for line in tinker_atoms:
		n = float(line[2])
		a = line[3]
		atom_dict[n] = a
	improper_atoms = tinker_impropers[:,0:4].astype('int')
	for q,(i,j,k,l) in enumerate(improper_atoms):
		fc,phi0 = tinker_impropers[q,-2:]
		fc = fc*4.184
		a0,a1,a2,a3 = atom_dict[i],atom_dict[j],atom_dict[k],atom_dict[l]
		s  = ' <Improper class1="%s" class2="%s" class3="%s" class4="%s"'%(a0,a1,a2,a3)
		s += ' phi0="%f" fc="%f"'%(phi0,fc)
		s += '/>'
		print(s)

def fix_wildcards(x):
	if x=="*":
		return ""
	else:
		return x

def load_vipar(filename):
	f = open(filename)
	s = f.readlines()
	sjoin = string.join(s)
	return eval(sjoin)

def lookup(atom, template):
	d = {}
	for k,[a,x,x,x] in enumerate(template["atoms"]):
		d[a] = k
	return d[atom]

angles = load_vipar(Root+"/angle_harm")
bonds = load_vipar(Root+"/stretch_harm")
impropers = load_vipar(Root+"/improper_harm")
torsions = load_vipar(Root+"/torsiontorsion_cmap")
ureybradley = load_vipar(Root+"/ureybradley_harm")
masses = load_vipar(Root+"/mass")
cmap = load_vipar(Root+"/cmap")
vdw = load_vipar(Root+"/vdw1")
vdw14 = load_vipar(Root+"/vdw1_14")
templates = load_vipar(Root+"/templates")
propers = load_vipar(Root+"/dihedral_trig")

external_names = ["$2","$4","$7"]

element_dict = {
"Na+":"Na",
"Mg2+":"Mg",
'Cs+':"Cs",
'Ca2+':"Ca",
"Cl-":"Cl",
"Zn2+":"Zn"
}

sigma_dict = {}
epsilon_dict = {}
sigma14_dict = {}
epsilon14_dict = {}

for vdwitem in vdw:
	name = vdwitem["type"][0]
	sigma = vdwitem["params"]["sigma"]
	epsilon = vdwitem["params"]["epsilon"]
	sigma_dict[name] = sigma
	epsilon_dict[name] = epsilon
	sigma14_dict[name] = sigma
	epsilon14_dict[name] = epsilon
	
for vdwitem in vdw14:
	name = vdwitem["type"][0]
	sigma = vdwitem["params"]["sigma"]
	epsilon = vdwitem["params"]["epsilon"]
	sigma14_dict[name] = sigma
	epsilon14_dict[name] = epsilon


print '<ForceField>'
mass_lookup = {}
new_element_dict = {}
for massitem in masses:
	name = massitem["type"][0]
	mass = massitem["params"]["amu"]
	try:
		element = element_dict[name]
	except KeyError:
		element = name[0:1]
	new_element_dict[name] = element
	mass_lookup[name] = mass

atom_properties = {}
vdw_properties = {}

print("<Residues>")
for name, template in templates.iteritems():
	print("""<Residue name ="%s">"""%name)
	for atom, element_number,charge,atom_class in template["atoms"]:
		atom_charge_class = "%s_%s"%(atom,name)
		print('<Atom name ="%s" type ="%s"/>'%(atom,atom_charge_class))
		atom_class = atom_class[0]
		sigma = sigma_dict[atom_class]
		epsilon = epsilon_dict[atom_class]
		sigma14 = sigma14_dict[atom_class]
		epsilon14 = epsilon14_dict[atom_class]
		atom_properties[atom_charge_class] = {"charge":charge,"sigma":sigma,"epsilon":epsilon,"sigma14":sigma14,"epsilon14":epsilon14,"class":atom_class}
		if vdw_properties.has_key(atom_class):
			if vdw_properties[atom_class] != {"sigma":sigma,"epsilon":epsilon,"sigma14":sigma14,"epsilon14":epsilon14}:
				raise Exception("Error: %s has incorrect VDW properties"%atom_class)
		vdw_properties[atom_class] = {"sigma":sigma,"epsilon":epsilon,"sigma14":sigma14,"epsilon14":epsilon14}
		
	if template.has_key("bonds"):
		for a0, a1 in template["bonds"]:
			if a0 and a1 not in external_names:
				b0 = lookup(a0,template)
				b1 = lookup(a1,template)
				print('<Bond from="%s" to="%s"/>'%(b0,b1))
		for a0, a1 in template["bonds"]:
			if a0 in external_names and a1 in external_names:
				print(a0,a1,"PASS")
			elif a0 in external_names:
				b = lookup(a1,template)
				print('<ExternalBond from="%s"/>'%(b))
			elif a1 in external_names:
				b = lookup(a0,template)
				print('<ExternalBond from="%s"/>'%(b))
	if template.has_key("cmap"):
		cmap_atoms = template["cmap"][1:-1]
		
	print("</Residue>")
print("</Residues>")

print ' <AtomTypes>'
for name,atomitem in atom_properties.iteritems():
	charge,sigma,epsilon,mass_class = atomitem["charge"],atomitem["sigma"],atomitem["epsilon"],atomitem["class"]
	element = new_element_dict[mass_class]
	mass = mass_lookup[mass_class]
	print '  <Type name="%s" class="%s" element="%s" mass="%s"/>' % (name,mass_class,element,mass)

print ' </AtomTypes>'



print ' <HarmonicBondForce>'
for bonditem in bonds:
	a0,a1 = bonditem["type"]
	r0 = bonditem["params"]["r0"] * length_conversion
	fc = bonditem["params"]["fc"] * energy_conversion
	fc *= 2.0 * length_conversion**-2.
	print '  <Bond class1="%s" class2="%s" length="%g" k="%f"/>' % (a0,a1,r0, fc)
print ' </HarmonicBondForce>'

print ' <HarmonicAngleForce>'
for angleitem in angles:
	a0,a1,a2 = angleitem["type"]
	theta0 = angleitem["params"]["theta0"] * angle_conversion
	fc = angleitem["params"]["fc"] * energy_conversion
	fc *= 2
	print '  <Angle class1="%s" class2="%s" class3="%s" angle="%f" k="%f"/>' % (a0,a1,a2,theta0, fc)
print ' </HarmonicAngleForce>'


torsion_dict = {}

print ' <PeriodicTorsionForce>'
for torsionitem in propers:
	a0,a1,a2,a3 = torsionitem["type"]
	b0,b1,b2,b3 = fix_wildcards(a0),fix_wildcards(a1),fix_wildcards(a2),fix_wildcards(a3)
	phi0 = torsionitem["params"]["phi0"] * np.pi/180.
	fclist = [torsionitem["params"]["fc%d"%i] * energy_conversion * 1.0 for i in range(7)]
	
	fclist[0] = (fclist[0] - sum(fclist[1:])) / (1+ np.cos(-1*phi0))
	
	if not torsion_dict.has_key((b0,b1,b2,b3)):
		torsion_dict[b0,b1,b2,b3] = {}
	
	torsion_dict[b0,b1,b2,b3][phi0] = fclist
	
for (b0,b1,b2,b3) in torsion_dict.keys():
	current_parms = []
	for theta,fclist in torsion_dict[b0,b1,b2,b3].iteritems():
		for k,f in enumerate(fclist):
			current_parms.append([f,k,theta])
	
	s  = ' <Proper class1="%s" class2="%s" class3="%s" class4="%s"'%(b0,b1,b2,b3)
	for i, (f,k,theta) in enumerate(current_parms):
		s += ' k%d="%f" periodicity%d="%d" phase%d="%f" '%(i+1,f,i+1,k,i+1,theta)
	
	s += "/>"

	print(s)
print ' </PeriodicTorsionForce>'


""" <CustomTorsionForce energy="fc*(theta - phi0)^2.">
<PerTorsionParameter name="phi0"/>
<PerTorsionParameter name="fc"/>
"""

#Use torsions that are derived from Tinker, rather than Desmond.
#This avoids double-counting some things
print "<ImproperHarmonicTorsionForce>"
ProcessTinkerImpropers()
print "</ImproperHarmonicTorsionForce>"
#print "</CustomTorsionForce>"

"""
for torsionitem in impropers:
	a0,a1,a2,a3 = torsionitem["type"]
	b0,b1,b2,b3 = fix_wildcards(a0),fix_wildcards(a1),fix_wildcards(a2),fix_wildcards(a3)
	phi0 = torsionitem["params"]["phi0"] * np.pi/180.
	fc = torsionitem["params"]["fc"] * energy_conversion

	#Fixes strange viparr bug IMHO
	switch_list = [("O","","","C"),("NC2","","","C"),("OC","","","CC")]#,("O","NH2","CT2","CC"),("O","CT2","NH2","CC")]
	if (b0,b1,b2,b3) in switch_list or (b0=="O" and b3=="CC"):
		b0,b1,b2,b3 = b3,b2,b1,b0
	
	s  = ' <Improper class1="%s" class2="%s" class3="%s" class4="%s"'%(b0,b1,b2,b3)

	s += ' phi0="%f" fc="%f"'%(phi0,fc)
	s += '/>'
	#print(s)
	#print("Replace me with tinker-derived torsions")
"""

print """ <UreyBradleyForce>"""
for angleitem in ureybradley:
	a0,a1,a2 = angleitem["type"]
	b0,b1,b2 = fix_wildcards(a0),fix_wildcards(a1),fix_wildcards(a2)
	r0 = angleitem["params"]["r0"] * length_conversion
	fc = angleitem["params"]["fc"] * energy_conversion
	fc *= 1.0*length_conversion**-2. #No factor of 2.0 because energy is same formula as desmond
	s  = ' <Urey class1="%s" class2="%s" class3="%s" r0="%f" fc="%f"/>'%(b0,b1,b2,r0, fc)
	print(s)
print """ </UreyBradleyForce>"""

print(""" <CharmmNonbondedForce coulomb14scale="1.0" lj14scale="1.0">""")
for name,atomitem in atom_properties.iteritems():
	charge,sigma,epsilon,atomclass = atomitem["charge"],atomitem["sigma"]*length_conversion,atomitem["epsilon"]*energy_conversion,atomitem["class"]
	sigma14,epsilon14 = atomitem["sigma14"]*length_conversion,atomitem["epsilon14"]*energy_conversion
	print("""<Atom type="%s" charge="%f" sigma="%f" epsilon="%f" sigma14="%f" epsilon14="%f"/>"""%(name,charge,sigma,epsilon,sigma14,epsilon14))

#print("""</NonbondedForce>""")
print("""</CharmmNonbondedForce>""")

print("""<CMAPTorsionForce>""")
for mapnum,c in enumerate(cmap):
	c = np.array(c)
	n = int(len(c)**0.5)
	phi,psi,c = c.transpose()
	ilist = (n*(phi%360)/360.).astype('int')
	jlist = (n*(psi%360)/360.).astype('int')
	energies = np.zeros((n**2))
	for k in range(n**2):
		i,j = ilist[k],jlist[k]
		energies[i+n*j] = c[k]
	
	energies *= energy_conversion
	
	print("<Map>")
	print(str(energies)[1:-1])
	print("</Map>")

print("""<Torsion map="2" class1="C" class2="N" class3="CP1" class4="C" class5="NH1"/>""")
print("""<Torsion map="3" class1="C" class2="N" class3="CP1" class4="C" class5="N"/>""")
print("""<Torsion map="4" class1="C" class2="NH1" class3="CT2" class4="C" class5="NH1"/>""")
print("""<Torsion map="5" class1="C" class2="NH1" class3="CT2" class4="C" class5="N"/>""")

print("""</CMAPTorsionForce>""")

print("""</ForceField>""")
