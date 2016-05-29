from __future__ import print_function
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as unit
import sys
from datetime import datetime
from optparse import OptionParser

def timeIntegration(context, steps, initialSteps):
    """Integrate a Context for a specified number of steps, then return how many seconds it took."""
    context.getIntegrator().step(initialSteps) # Make sure everything is fully initialized
    context.getState(getEnergy=True)
    start = datetime.now()
    context.getIntegrator().step(steps)
    context.getState(getEnergy=True)
    end = datetime.now()
    elapsed = end -start
    return elapsed.seconds + elapsed.microseconds*1e-6

def runOneTest(testName, options):
    """Perform a single benchmarking simulation."""
    explicit = (testName in ('rf', 'pme', 'amoebapme'))
    amoeba = (testName in ('amoebagk', 'amoebapme'))
    hydrogenMass = None
    print()
    if amoeba:
        print('Test: %s (epsilon=%g)' % (testName, options.epsilon))
    elif testName == 'pme':
        print('Test: pme (cutoff=%g)' % options.cutoff)
    else:
        print('Test: %s' % testName)
    platform = mm.Platform.getPlatformByName(options.platform)
    
    # Create the System.
    
    if amoeba:
        constraints = None
        epsilon = float(options.epsilon)
        if explicit:
            ff = app.ForceField('amoeba2009.xml')
            pdb = app.PDBFile('5dfr_solv-cube_equil.pdb')
            cutoff = 0.7*unit.nanometers
            vdwCutoff = 0.9*unit.nanometers
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, vdwCutoff=vdwCutoff, constraints=constraints, ewaldErrorTolerance=0.00075, mutualInducedTargetEpsilon=epsilon, polarization=options.polarization)
        else:
            ff = app.ForceField('amoeba2009.xml', 'amoeba2009_gk.xml')
            pdb = app.PDBFile('5dfr_minimized.pdb')
            cutoff = 2.0*unit.nanometers
            vdwCutoff = 1.2*unit.nanometers
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=constraints, mutualInducedTargetEpsilon=epsilon, polarization=options.polarization)
        for f in system.getForces():
            if isinstance(f, mm.AmoebaMultipoleForce) or isinstance(f, mm.AmoebaVdwForce) or isinstance(f, mm.AmoebaGeneralizedKirkwoodForce) or isinstance(f, mm.AmoebaWcaDispersionForce):
                f.setForceGroup(1)
        dt = 0.002*unit.picoseconds
        integ = mm.MTSIntegrator(dt, [(0,2), (1,1)])
    else:
        if explicit:
            ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
            pdb = app.PDBFile('5dfr_solv-cube_equil.pdb')
            if testName == 'pme':
                method = app.PME
                cutoff = options.cutoff
            else:
                method = app.CutoffPeriodic
                cutoff = 1*unit.nanometers
        else:
            ff = app.ForceField('amber99sb.xml', 'amber99_obc.xml')
            pdb = app.PDBFile('5dfr_minimized.pdb')
            method = app.CutoffNonPeriodic
            cutoff = 2*unit.nanometers
        if options.heavy:
            dt = 0.005*unit.picoseconds
            constraints = app.AllBonds
            hydrogenMass = 4*unit.amu
        else:
            dt = 0.002*unit.picoseconds
            constraints = app.HBonds
            hydrogenMass = None
        system = ff.createSystem(pdb.topology, nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints, hydrogenMass=hydrogenMass)
        integ = mm.LangevinIntegrator(300*unit.kelvin, 91*(1/unit.picoseconds), dt)
    print('Step Size: %g fs' % dt.value_in_unit(unit.femtoseconds))
    properties = {}
    initialSteps = 5
    if options.device is not None and platform.getName() in ('CUDA', 'OpenCL'):
        properties['DeviceIndex'] = options.device
        if ',' in options.device or ' ' in options.device:
            initialSteps = 250
    if options.precision is not None and platform.getName() in ('CUDA', 'OpenCL'):
        properties['Precision'] = options.precision
    
    # Run the simulation.
    
    integ.setConstraintTolerance(1e-5)
    if len(properties) > 0:
        context = mm.Context(system, integ, platform, properties)
    else:
        context = mm.Context(system, integ, platform)
    context.setPositions(pdb.positions)
    context.setVelocitiesToTemperature(300*unit.kelvin)
    steps = 20
    while True:
        time = timeIntegration(context, steps, initialSteps)
        if time >= 0.5*options.seconds:
            break
        if time < 0.5:
            steps = int(steps*1.0/time) # Integrate enough steps to get a reasonable estimate for how many we'll need.
        else:
            steps = int(steps*options.seconds/time)
    print('Integrated %d steps in %g seconds' % (steps, time))
    print('%g ns/day' % (dt*steps*86400/time).value_in_unit(unit.nanoseconds))

# Parse the command line options.

parser = OptionParser()
platformNames = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
parser.add_option('--platform', dest='platform', choices=platformNames, help='name of the platform to benchmark')
parser.add_option('--test', dest='test', choices=('gbsa', 'rf', 'pme', 'amoebagk', 'amoebapme'), help='the test to perform: gbsa, rf, pme, amoebagk, or amoebapme [default: all]')
parser.add_option('--pme-cutoff', default='0.9', dest='cutoff', type='float', help='direct space cutoff for PME in nm [default: 0.9]')
parser.add_option('--seconds', default='60', dest='seconds', type='float', help='target simulation length in seconds [default: 60]')
parser.add_option('--polarization', default='mutual', dest='polarization', choices=('direct', 'extrapolated', 'mutual'), help='the polarization method for AMOEBA: direct, extrapolated, or mutual [default: mutual]')
parser.add_option('--mutual-epsilon', default='1e-5', dest='epsilon', type='float', help='mutual induced epsilon for AMOEBA [default: 1e-5]')
parser.add_option('--heavy-hydrogens', action='store_true', default=False, dest='heavy', help='repartition mass to allow a larger time step')
parser.add_option('--device', default=None, dest='device', help='device index for CUDA or OpenCL')
parser.add_option('--precision', default='single', dest='precision', choices=('single', 'mixed', 'double'), help='precision mode for CUDA or OpenCL: single, mixed, or double [default: single]')
(options, args) = parser.parse_args()
if len(args) > 0:
    parser.error('Unknown argument: '+args[0])
if options.platform is None:
    parser.error('No platform specified')
print('Platform:', options.platform)
if options.platform in ('CUDA', 'OpenCL'):
    print('Precision:', options.precision)
    if options.device is not None:
        print('Device:', options.device)

# Run the simulations.

if options.test is None:
    for test in ('gbsa', 'rf', 'pme', 'amoebagk', 'amoebapme'):
        try:
            runOneTest(test, options)
        except Exception as ex:
            print('Test failed: %s' % ex.message)
else:
    runOneTest(options.test, options)
