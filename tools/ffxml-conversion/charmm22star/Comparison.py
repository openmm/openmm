from simtk import unit as unit
import numpy as np
import DMSReader
import OpenMM_Data

eps = 5E-5

class Comparison():
	def __init__(self,dms_filename,pdb,system):
		self.omm = OpenMM_Data.OpenMM_Data(system)
		self.dms = DMSReader.DMS_Data(dms_filename)
		
		self.atoms = np.array([a.name for a in pdb.topology.atoms()])
		self.na = len(self.atoms)
	
	def old_construct_vdw_dict(self):
		self.sigma_dict_mae = {}
		self.epsilon_dict_mae = {}
		
		for i, vec in enumerate(self.vdw_mae):
			key = vec[1]
			sigma = float(vec[3])
			epsilon = float(vec[4])
			self.sigma_dict_mae[key] = sigma
			self.epsilon_dict_mae[key] = epsilon
			
	def compare(self):
		self.compare_bonds_and_urey()
		self.compare_angles()
		self.compare_dihedral_trigs()
		self.compare_improper_harms()
		self.compare_nonbonded()
		self.compare_exceptions()
		self.compare_cmap()

	def compare_bonds_and_urey(self):
		for key,(r0,fc) in self.dms.bonds.iteritems():
			if key not in self.omm.bonds.keys():
				print(key,r0,fc)
				raise(Exception("Bond / Urey Pattern Different"))

		for key,(r0,fc) in self.omm.bonds.iteritems():
			if key not in self.dms.bonds.keys():
				print(key,r0,fc)
				raise(Exception("Bond / Urey Pattern Different"))

		for key,(r0,fc) in self.omm.bonds.iteritems():
			r0_2,fc_2 = self.dms.bonds[key]
			if abs(r0 - r0_2) > eps or abs(fc - fc_2) > eps:
				print(key,r0,fc,r0_2,fc_2)
				raise(Exception("Bond / Urey Term Different"))

	def compare_angles(self):
		for key,(r0,fc) in self.dms.angles.iteritems():
			if key not in self.omm.angles.keys():
				print(key,r0,fc)
				raise(Exception("Angle Pattern Different"))

		for key,(r0,fc) in self.omm.angles.iteritems():
			if key not in self.dms.angles.keys():
				print(key,r0,fc)
				raise(Exception("Angle Pattern Different"))

		for key,(r0,fc) in self.omm.angles.iteritems():
			r0_2,fc_2 = self.dms.angles[key]
			if abs(r0 - r0_2) > eps or abs(fc - fc_2) > eps:
				print(key,r0,fc,r0_2,fc_2)
				raise(Exception("Angle Term Different"))


	def compare_dihedral_trigs(self):
		for key,val in self.dms.dihedral_trigs.iteritems():
			if key not in self.omm.dihedral_trigs.keys():
				print(key)
				raise(Exception("dihedral_trig Pattern Different"))

		for key,val in self.omm.dihedral_trigs.iteritems():
			if key not in self.dms.dihedral_trigs.keys():
				print(key)
				raise(Exception("dihedral_trig Pattern Different"))

		for (p0,p1,p2,p3),val in self.omm.dihedral_trigs.iteritems():
			for phi0,fc_omm in val.iteritems():
				fc_dms = self.dms.dihedral_trigs[p0,p1,p2,p3][phi0]
				tmp = fc_omm[:]*1.0
				tmp[0] = tmp.sum() + tmp[0]*np.cos(phi0)
				if np.linalg.norm(fc_dms - tmp) > eps:
					print(p0,p1,p2,p3)
					print(fc_dms)
					print(tmp)
					raise(Exception("dihedral_trig Term Different"))

		for (p0,p1,p2,p3),val in self.dms.dihedral_trigs.iteritems():
			for phi0,fc_dms in val.iteritems():
				fc_omm = self.omm.dihedral_trigs[p0,p1,p2,p3][phi0]
				tmp = fc_omm[:]*1.0
				tmp[0] = tmp.sum() + tmp[0]*np.cos(phi0)
				if np.linalg.norm(fc_dms - tmp) > eps:
					print(p0,p1,p2,p3)
					print(fc_dms)
					print(tmp)
					raise(Exception("dihedral_trig Term Different"))		
		

	def compare_improper_harms(self):
		for (p0,p1,p2,p3),(phi0,fc) in self.dms.improper_harms.iteritems():
			phi0_mm,fc_mm = self.omm.improper_harms[p0,p1,p2,p3]
			if abs(phi0-phi0_mm) > eps or abs(fc-fc_mm) > eps:
					print(p0,p1,p2,p3,phi0,fc,phi0_mm,fc_mm)
					raise(Exception("improper_harm Term Different"))

	def compare_nonbonded(self):
		for fid,q in self.dms.charges.iteritems():
			q_omm = self.omm.charges[fid]
			if abs(q_omm - q) > eps:
				print(fid,q,qomm)
				raise(Exception("Charges Different"))

		for fid,q in self.dms.sigmas.iteritems():
			q_omm = self.omm.sigmas[fid]
			if abs(q_omm - q) > eps:
				print(fid,q,qomm)
				raise(Exception("Sigmas Different"))

		for fid,q in self.dms.epsilons.iteritems():
			q_omm = self.omm.epsilons[fid]
			if abs(q_omm - q) > eps:
				print(fid,q,q_omm)
				raise(Exception("Epsilons Different"))

	def compare_exceptions(self):
		for (i,j), (sigma,epsilon,q) in self.dms.exceptions.iteritems():
			s2,e2,q2 = self.omm.exceptions[i,j]
			if abs(s2 - sigma) > eps or abs(e2 - epsilon) > eps or abs(q2 - q) > eps:
				print(i,j,sigma,epsilon,q,s2,e2,q2)
				raise(Exception("Exceptions Different"))

		for (i,j), (sigma,epsilon,q) in self.omm.exceptions.iteritems():
			
			if self.dms.exceptions.has_key((i,j)):
				s2,e2,q2 = self.dms.exceptions[i,j]
			elif q == 0 and epsilon == 0:
				continue
			else:
				s2i = self.dms.sigmas[i]
				s2j = self.dms.sigmas[j]
				e2i = self.dms.epsilons[i]
				e2j = self.dms.epsilons[j]
				q2i = self.dms.charges[i]
				q2j = self.dms.charges[j]
				s2 = 0.5*(s2i + s2j)
				e2 = (e2i*e2j)**0.5
				q2 = q2i*q2j
			
			if abs(s2 - sigma) > eps or abs(e2 - epsilon) > eps or abs(q2 - q) > eps:
				print(i,j,sigma,epsilon,q,s2,e2,q2)
				raise(Exception("Exceptions Different"))

	def compare_cmap(self):
		pass
