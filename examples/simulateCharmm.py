from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Read the PSF
psf = CharmmPsfFile('ala_ala_ala.psf')

# Get the coordinates from the PDB
pdb = PDBFile('ala_ala_ala.pdb')

# Read and condense the parameter set
params = CharmmParameterSet('charmm22.rtf', 'charmm22.par')

# Instantiate the system
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
                          nonbondedCutoff=None)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.getPositions())
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('ch_output.pdb', 1000))
simulation.reporters.append(
        StateDataReporter(stdout, 1000, step=True, potentialEnergy=True,
                          temperature=True)
)
simulation.step(10000)
