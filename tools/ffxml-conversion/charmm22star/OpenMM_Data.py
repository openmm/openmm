import simtk.openmm
from simtk import unit as unit
import numpy as np
	
class OpenMM_Data():
	def __init__(self,system):
		self.prepare_openmm(system)

	def prepare_openmm(self,system):
		num_forces = system.getNumForces()
		forces = [system.getForce(i) for i in range(num_forces)]
		for i in range(num_forces):
			if forces[i].__class__ == simtk.openmm.openmm.HarmonicBondForce:
				self.f_bond = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.HarmonicAngleForce:
				self.f_angle = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.PeriodicTorsionForce:
				self.f_torsion = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.RBTorsionForce:
				self.f_RB = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.NonbondedForce:
				self.f_NB = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.CustomCompoundBondForce:
				self.f_urey = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.CustomTorsionForce:
				self.f_improper_harm = forces[i]
			elif forces[i].__class__ == simtk.openmm.openmm.CMAPTorsionForce:
				self.f_cmap = forces[i]

		self.get_bonds(self.f_bond)
		self.get_angles(self.f_angle)
		self.get_dihedral_trigs(self.f_torsion)
		self.get_ureys(self.f_urey)
		self.get_improper_harms(self.f_improper_harm)
		self.get_nonbonded(self.f_NB)						
		self.get_cmap(self.f_cmap)
	def get_bonds(self,F):
		self.bonds = {}
		for m in range(F.getNumBonds()):
			i,j,r0,fc = F.getBondParameters(m)
			r0 = r0.value_in_unit(unit.angstrom)
			fc = fc.value_in_unit(unit.kilocalories_per_mole / unit.angstrom**2)
			self.bonds[i,j] = (r0,fc)
			self.bonds[j,i] = (r0,fc)

	def get_angles(self,F):
		self.angles = {}
		for m in range(F.getNumAngles()):
			i,j,k,phi0,fc = F.getAngleParameters(m)
			phi0 = phi0.value_in_unit(unit.degrees)
			fc = fc.value_in_unit(unit.kilocalories_per_mole / (unit.radian**2))
			self.angles[i,j,k] = (phi0,fc)
			self.angles[k,j,i] = (phi0,fc)

	def get_improper_harms(self,F):
		self.improper_harms = {}
		for m in range(F.getNumTorsions()):
			i,j,k,l,(phi0,fc) = F.getTorsionParameters(m)
			fc = fc / 4.184
			s = tuple(np.unique([i,j,k,l]))
			if fc > 0:
				self.improper_harms[s] = (phi0,fc)
		
	def get_ureys(self,F):
		"""Note: we add urey-bradleys to the bond dictionary because they both share the same functional form"""
		for m in range(F.getNumBonds()):
			(i,j,k),(r0,fc) = F.getBondParameters(m)
			r0 *= 10.
			fc /= (4.184*100*0.5)
			self.bonds[i,k] = (r0,fc)
			self.bonds[k,i] = (r0,fc)

	def get_dihedral_trigs(self,F):
		self.dihedral_trigs = {}
		for m in range(F.getNumTorsions()):
			i,j,k,l,periodicity,phase,strength = F.getTorsionParameters(m)
			s = strength.value_in_unit(unit.kilocalories_per_mole)
			theta = np.float16(phase.value_in_unit(unit.degree))
			
			if not self.dihedral_trigs.has_key((i,j,k,l)):
				self.dihedral_trigs[(i,j,k,l)] = {}
				self.dihedral_trigs[(l,k,j,i)] = {}
				
			if not self.dihedral_trigs[(i,j,k,l)].has_key(theta):
				z = np.zeros(7)
				self.dihedral_trigs[(i,j,k,l)][theta] = z
				self.dihedral_trigs[(l,k,j,i)][theta] = z
				z[periodicity] = s
			else:
				self.dihedral_trigs[(i,j,k,l)][theta][periodicity] = s
				self.dihedral_trigs[(l,k,j,i)][theta][periodicity] = s
	
	def get_nonbonded(self,F):
		self.charges = {}
		self.sigmas = {}
		self.epsilons = {}
		
		for m in range(F.getNumParticles()):
			q,s,e = F.getParticleParameters(m)
			q = q.value_in_unit(unit.elementary_charge)
			s = s.value_in_unit(unit.angstrom)
			e = e.value_in_unit(unit.kilocalories_per_mole)
			self.charges[m] = q
			self.sigmas[m] = s
			self.epsilons[m] = e
		
		self.exceptions = {}
		for m in range(F.getNumExceptions()):
			i,j,q,s,e = F.getExceptionParameters(m)
			q = q.value_in_unit(unit.elementary_charge*unit.elementary_charge)
			s = s.value_in_unit(unit.angstrom)
			e = e.value_in_unit(unit.kilocalories_per_mole)
			self.exceptions[(i,j)] = (s,e,q)
			self.exceptions[(j,i)] = (s,e,q)

	def get_cmap(self,F):
		self.cmap_ids = {}
		
		for m in range(F.getNumTorsions()):
			fid,q0,q1,q2,q3,q4,q5,q6,q7 = F.getTorsionParameters(m)
			self.cmap_ids[q0,q1,q2,q3,q4,q5,q6,q7] = fid
		
		self.cmap_grids = {}
		for m in range(F.getNumMaps()):
			ngrid,grid = F.getMapParameters(m)
			self.cmap_grids[m] = grid
