import numpy as np
import sqlite3

eps = 1E-5

class DMS_Data():

	def __init__(self,filename):
		self.conn = sqlite3.connect(filename)
		self.c = self.conn.cursor()
		self.get_nonbonded()
		self.get_bonds()
		self.get_angles()
		self.get_dihedral_trigs()
		self.get_improper_harms()
		self.get_exceptions()
		self.get_cmap()

	def get_exceptions(self):
		self.c.execute("SELECT * FROM pair_12_6_es_param")
		ans = self.c.fetchall()
		
		self.exception_types = {}
		for row in ans:
			aij,bij,qij,fid,typekey,comment,ff = row
			self.exception_types[fid] = (aij,bij,qij)

		self.c.execute("SELECT * FROM pair_12_6_es_term")
		ans = self.c.fetchall()
		
		self.exceptions = {}
		for row in ans:
			p0,p1,fid,moiety = row
			a,b,q = self.exception_types[fid]
			sigma = (a/b)**(1.0/6.0)
			epsilon = b*b / (4.0*a)
			self.exceptions[p0,p1] = (sigma,epsilon,q)
			self.exceptions[p1,p0] = (sigma,epsilon,q)
			
	def get_nonbonded(self):
		
		self.c.execute("SELECT * FROM nonbonded_param")
		ans = self.c.fetchall()
		self.nonbonded_types = {}

		for row in ans:
			sigma,epsilon,fid,typekey,comment,ff = row
			self.nonbonded_types[fid] = (sigma,epsilon)
		
		self.c.execute("SELECT * FROM particle")
		ans = self.c.fetchall()

		self.resnames = {}
		self.names = {}
		self.charges = {}
		self.resids = {}
		self.sigmas = {}
		self.epsilons = {}
		
		for row in ans:
			key, anum, chain, resname, resid, name, x,y,z,vx,vy,vz,mass,charge,nbtype,grp_temperature,grp_energy,grp_ligand,grp_bias,grp_frozen,occupancy,bfactor,formal_charge = row
			self.resnames[key] = resname
			self.resids[key] = resid
			self.charges[key] = charge
			self.names[key] = name
			sigma,epsilon = self.nonbonded_types[nbtype]
			self.sigmas[key] = sigma
			self.epsilons[key] = epsilon

	def get_bonds(self):
		self.c.execute("SELECT * FROM stretch_harm_term")
		ans = self.c.fetchall()
		self.bond_param_ids = {}
		for row in ans:
			p0,p1,param,constrained,moiety = row
			self.bond_param_ids[(p0,p1)] = param
		
		self.c.execute("SELECT * FROM stretch_harm_param")
		ans = self.c.fetchall()
		self.bond_params = {}
		for row in ans:
			r0,fc,fid,typekey,comment,ff = row
			fc = fc * 2.0
			self.bond_params[fid] = (r0,fc)
		
		self.bonds = {}
		for key in self.bond_param_ids.keys():
			self.bonds[key] = self.bond_params[self.bond_param_ids[key]]
			p0,p1 = key
			key1 = (p1,p0)
			self.bonds[key1] = self.bond_params[self.bond_param_ids[key]]

	def get_angles(self):
		self.c.execute("SELECT * FROM angle_harm_term")
		ans = self.c.fetchall()
		self.angle_param_ids = {}
		for row in ans:
			p0,p1,p2,param,constrained,moiety = row
			self.angle_param_ids[(p0,p1,p2)] = param
		
		self.c.execute("SELECT * FROM angle_harm_param")
		ans = self.c.fetchall()
		self.angle_params = {}
		for row in ans:
			r0,fc,fid,typekey,comment,ff = row
			fc = fc * 2.0
			self.angle_params[fid] = (r0,fc)
		
		self.angles = {}
		for key in self.angle_param_ids.keys():
			self.angles[key] = self.angle_params[self.angle_param_ids[key]]
			p0,p1,p2 = key
			key1 = (p2,p1,p0)
			self.angles[key1] = self.angle_params[self.angle_param_ids[key]]


	def get_dihedral_trigs(self):
		self.c.execute("SELECT * FROM dihedral_trig_term")
		ans = self.c.fetchall()
		self.dihedral_trig_param_ids = {}
		for row in ans:
			p0,p1,p2,p3,param,moiety = row
			key = (p0,p1,p2,p3)
			if self.dihedral_trig_param_ids.has_key(key):
				self.dihedral_trig_param_ids[key].append(param)
			else:
				self.dihedral_trig_param_ids[key] = [param]
		
		self.c.execute("SELECT * FROM dihedral_trig_param")
		ans = self.c.fetchall()
		self.dihedral_trig_params = {}
		for row in ans:
			phi0,fc0,fc1,fc2,fc3,fc4,fc5,fc6,fid,typekey,comment,ff = row
			self.dihedral_trig_params[fid] = (phi0,fc0,fc1,fc2,fc3,fc4,fc5,fc6)
		
		self.dihedral_trigs = {}
		for key in self.dihedral_trig_param_ids.keys():
			p0,p1,p2,p3 = key
			key2 = p3,p2,p1,p0
			for term in self.dihedral_trig_param_ids[key]:
				(phi0,fc0,fc1,fc2,fc3,fc4,fc5,fc6) = self.dihedral_trig_params[term]
				parms = np.array([fc0,fc1,fc2,fc3,fc4,fc5,fc6])
				if np.linalg.norm(parms) > eps:
					if not self.dihedral_trigs.has_key(key):
						self.dihedral_trigs[key] = {}
						self.dihedral_trigs[key2] = {}
					self.dihedral_trigs[key][phi0] = parms
					self.dihedral_trigs[key2][phi0] = parms
		
	def get_improper_harms(self):
		self.c.execute("SELECT * FROM improper_harm_term")
		ans = self.c.fetchall()
		self.improper_harm_param_ids = {}
		for row in ans:
			p0,p1,p2,p3,param,moiety = row
			key = np.unique([p0,p1,p2,p3])
			key = tuple(key)
			if not self.improper_harm_param_ids.has_key(key):
				self.improper_harm_param_ids[key] = []
			
			self.improper_harm_param_ids[key].append(param)
		
		self.c.execute("SELECT * FROM improper_harm_param")
		ans = self.c.fetchall()
		self.improper_harm_params = {}
		for row in ans:
			phi0,fc,fid,typekey,comment,ff = row
			self.improper_harm_params[fid] = (phi0,fc)
		
		self.improper_harms = {}
		for key,fidlist in self.improper_harm_param_ids.iteritems():
			philist = []
			for fid in fidlist:
				phi0,fc = self.improper_harm_params[fid]
				if fc > 0.:
					if not self.improper_harms.has_key(key):
						self.improper_harms[key] = (phi0,0.0)
					oldfc = self.improper_harms[key][1]
					self.improper_harms[key] = (phi0,fc + oldfc)
					philist.append(phi0)
				if np.std(philist) > eps:
					raise(Exception("Putting multiple improper harmonics on single atom set"))
		
	def get_cmap(self):
		self.c.execute("SELECT * FROM torsiontorsion_cmap_term")
		ans = self.c.fetchall()
		self.cmap_param_ids = {}
		for row in ans:
			p0,p1,p2,p3,p4,p5,p6,p7,param,moiety = row
			key = p0,p1,p2,p3,p4,p5,p6,p7
			if not self.cmap_param_ids.has_key(key):
				self.cmap_param_ids[key] = []
			
			self.cmap_param_ids[key].append(param)

		self.c.execute("SELECT * FROM torsiontorsion_cmap_param")
		ans = self.c.fetchall()
		self.cmap_param_string = {}
		for row in ans:
			cmapid,fid,typekey,comment,ff = row
			self.cmap_param_string[fid] = (cmapid)

		self.cmap_grids = {}
		for key,cmapstring in self.cmap_param_string.iteritems():
			
			self.c.execute("SELECT * FROM %s"%cmapstring)
			ans = self.c.fetchall()
			philist = []
			psilist = []
			enlist = []
			for row in ans:
				phi,psi,energy = row
				philist.append(phi)
				psilist.append(psi)
				enlist.append(energy)
			self.cmap_grids[key] = (np.array(philist),np.array(psilist),np.array(enlist))
