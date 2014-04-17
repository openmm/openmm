from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr

# Define a user-interface
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-t', '--top', dest='top', default=None,
                  help='CHARMM RTF file to use in the simulation.')
parser.add_option('-p', '--param', dest='par', default=None,
                  help='CHARMM parameter file to use in the simulation.')
opt, arg = parser.parse_args()

if arg:
    stderr.write('Unexpected arguments: %s' % ', '.join(arg) + '\n')
    exit(parser.print_help() or 1)
if opt.top is None or opt.par is None:
    stderr.write('You must provide a top AND parameter file\n')
    exit(parser.print_help() or 1)

# Read the PSF
psf = CharmmPsfFile('ala_ala_ala.psf')

# Get the coordinates from the PDB
pdb = PDBFile('ala_ala_ala.pdb')

# Read and condense the parameter set
params = CharmmParameterSet(opt.top, opt.par)

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
