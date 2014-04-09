from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Read the PSF
psf = CharmmPSF.load_from_psf('ala_ala_ala.psf')

# Get the coordinates from the PDB
pdb = PDBFile('ala_ala_ala.pdb')

# Read and condense the parameter set
params = CharmmParameterSet.load_set(tfile='charmm22.rtf',
                                     pfile='charmm22.par').condense()

# Load the parameter set !! do not forget this !!
psf.load_parameters(params)

# Instantiate the system
system = psf.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=None)
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
