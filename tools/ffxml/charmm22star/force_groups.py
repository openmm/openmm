import string
import numpy as np

def create_force_groups(simulation):
	system = simulation.context.getSystem()
	n  = system.getNumForces()
	forces = []
	for i in range(n):
		f = system.getForce(i)
		f.setForceGroup(i)
		forces.append(f)
	
	simulation.context.reinitialize()
	return forces
	
def get_energies(simulation):
	context = simulation.context
	system = simulation.context.getSystem()
	n  = system.getNumForces()
	bits = np.zeros(n,'int')
	energies = []
	
	state = context.getState(getEnergy=True)
	total_energy = state.getPotentialEnergy()
	
	for i in range(n):
		"""groups = bits.copy()
		groups[i] = 1.
		groups = [str(g) for g in groups]
		groups = "0b"+string.join(groups,"")
		print(i,groups)
		groups = int(groups,2)
		"""
		groups = (1<<i)
		state = context.getState(getEnergy=True,groups=groups)
		energy = state.getPotentialEnergy()
		energies.append(energy)
	
	return total_energy,energies
