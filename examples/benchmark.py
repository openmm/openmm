from __future__ import print_function
import openmm.app as app
import openmm as mm
import openmm.unit as unit
from datetime import datetime, timezone
import argparse
import os

def cpuinfo():
    """Return CPU info"""
    import platform, subprocess, re
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        command ="sysctl -n machdep.cpu.brand_string"
        return subprocess.check_output(command, shell=True, text=True).strip()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True, text=True).strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub( ".*model name.*: ", "", line,1)
    return ""

def appendTestResult(filename=None, test_result=None, system_info=None):
    """Append a test result to a JSON or YAML file.

    TODO: The efficiency of this could be improved.

    Parameters
    ----------
    filename : str, optional, default=None
        The filename to append a result to, ending either in .yaml or .json
        If None, no result is written
    test_result : dict, optional, default=None
        The test result to append to the 'benchmarks' blcok
    system_info : dict, optional, default=None
        System info to append to the 'system' block
    """
    # Do nothing if filename is None
    if filename is None:
        return

    all_results = { 'benchmarks' : list() }
    if system_info is not None:
        all_results['system'] = system_info

    if filename.endswith('.yaml'):
        # Append test result to a YAML file, creating if file does not exist.
        import yaml
        if os.path.exists(filename):
            with open(filename, 'rt') as infile:
                all_results = yaml.safe_load(infile)
        if test_result is not None:
            all_results['benchmarks'].append(test_result)
        with open(filename, 'wt') as outfile:
            yaml.dump(all_results, outfile, sort_keys=False)
    elif filename.endswith('.json'):
        # Append test result to a JSON file, creating if file does not exist.
        import json
        if os.path.exists(filename):
            with open(filename, 'rt') as infile:
                all_results = json.load(infile)
        if test_result is not None:
            all_results['benchmarks'].append(test_result)
        with open(filename, 'wt') as outfile:
            json.dump(all_results, outfile, sort_keys=False, indent=4)
    else:
        raise ValueError('--output filename must end with .json or .yaml')

def printTestResult(test_result, options):
    """Render a test result

    Parameters
    ----------
    test_result : dict
        The test result
    options :
        Options structure
    """
    if options.style == 'simple':
        for (key, value) in test_result.items():
            print(f'{key}: {value}')
        print('')
    elif options.style == 'table':
        print('%-18s%-12s%-14s%-15s%-10g%-11s%-11s%-g' %
              (test_result['test'],
               test_result['precision'],
               test_result['constraints'],
               test_result['hydrogen_mass'],
               test_result['timestep_in_fs'],
               test_result['ensemble'],
               test_result['platform'],
               test_result['ns_per_day']))
    else:
        raise ValueError(f"style '{style}' must be one of ['simple', 'table']")

def timeIntegration(context, steps, initialSteps):
    """Integrate a Context for a specified number of steps, then return how many seconds it took."""
    context.getIntegrator().step(initialSteps) # Make sure everything is fully initialized
    context.getState(energy=True)
    start = datetime.now()
    context.getIntegrator().step(steps)
    context.getState(energy=True)
    end = datetime.now()
    elapsed = end-start
    return elapsed.seconds + elapsed.microseconds*1e-6

def downloadAmberSuite():
    """Download and extract Amber benchmark to Amber20_Benchmark_Suite/ in current directory."""
    dirname = 'Amber20_Benchmark_Suite'
    url = 'https://ambermd.org/Amber20_Benchmark_Suite.tar.gz'
    if not os.path.exists(dirname):
        import urllib.request
        print('Downloading', url)
        filename, headers = urllib.request.urlretrieve(url, filename='Amber20_Benchmark_Suite.tar.gz')
        import tarfile
        print('Extracting', filename)
        tarfh = tarfile.open(filename, 'r:gz')
        tarfh.extractall(path=dirname)
    return dirname

import functools
@functools.lru_cache(maxsize=None)
def retrieveTestSystem(testName, pme_cutoff=0.9, bond_constraints='hbonds', polarization='mutual', epsilon=1e-5):
    """Retrieve a benchmark system

    Parameters
    ----------
    testName : str
        The name of the test
    pme_cutoff : float
        PME cutoff, in nm
    bond_constraints : str, optional, default='hbonds'
        'hbonds' : use app.HBonds and set H mass to 1.5*amu
        'allbonds' : use app.AllBonds and set H mass to 4*amu
    polarization : str, optional, default='mutual'
        Polarization scheme for Amoeba
    epsilon : str or float, optional, default=1e-5
        mutualInducedTargetEpsilon for Amoeba

    Returns
    -------
    system : openmm.System
        The test system object
    positions : openmm.unit.Quantity with shape (natoms,3)
        The initial positions
    test_parameters : dict of str : str
        Special test parameters to report in test results

    """
    explicit = (testName not in ('gbsa', 'amoebagk'))
    amoeba = (testName in ('amoebagk', 'amoebapme'))
    apoa1 = testName.startswith('apoa1')
    amber = (testName.startswith('amber'))
    hydrogenMass = None

    # Create dictionary of test parameters
    test_parameters = dict()

    # Store the test name
    test_parameters['test'] = testName

    # Create the System.
    if amoeba:
        constraints = None
        test_parameters['epsilon'] = epsilon
        epsilon = float(epsilon)
        if explicit:
            ff = app.ForceField('amoeba2009.xml')
            pdb = app.PDBFile('5dfr_solv-cube_equil.pdb')
            cutoff = 0.7*unit.nanometers
            vdwCutoff = 0.9*unit.nanometers
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, vdwCutoff=vdwCutoff, constraints=constraints, ewaldErrorTolerance=0.00075, mutualInducedTargetEpsilon=epsilon, polarization=polarization)
        else:
            ff = app.ForceField('amoeba2009.xml', 'amoeba2009_gk.xml')
            pdb = app.PDBFile('5dfr_minimized.pdb')
            system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=constraints, mutualInducedTargetEpsilon=epsilon, polarization=polarization)
        for f in system.getForces():
            if isinstance(f, mm.AmoebaMultipoleForce) or isinstance(f, mm.AmoebaVdwForce) or isinstance(f, mm.AmoebaGeneralizedKirkwoodForce) or isinstance(f, mm.AmoebaWcaDispersionForce):
                f.setForceGroup(1)
        positions = pdb.positions
        test_parameters['constraints'] = 'None'
        test_parameters['hydrogen_mass'] = '1'
    elif amber:
        dirname = downloadAmberSuite()
        names = {'amber20-dhfr':'JAC', 'amber20-cellulose':'Cellulose', 'amber20-stmv':'STMV'}
        fileName = names[testName]
        prmtop = app.AmberPrmtopFile(os.path.join(dirname, f'PME/Topologies/{fileName}.prmtop'))
        inpcrd = app.AmberInpcrdFile(os.path.join(dirname, f'PME/Coordinates/{fileName}.inpcrd'))
        positions = inpcrd.positions
        method = app.PME
        constraints = app.HBonds
        cutoff = pme_cutoff
        system = prmtop.createSystem(nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints)
        if inpcrd.boxVectors is not None:
            system.setDefaultPeriodicBoxVectors(*inpcrd.boxVectors)
        test_parameters['cutoff'] = cutoff
        test_parameters['constraints'] = 'HBonds'
        test_parameters['hydrogen_mass'] = '1.5'
    else:
        if apoa1:
            ff = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/lipid17.xml', 'amber14/tip3p.xml')
            pdb = app.PDBFile('apoa1.pdb')
            if testName == 'apoa1pme':
                method = app.PME
                cutoff = pme_cutoff
            elif testName == 'apoa1ljpme':
                method = app.LJPME
                cutoff = pme_cutoff
            else:
                # Reaction field uses hard-coded cutoff
                method = app.CutoffPeriodic
                cutoff = 1.0 # nanometers ; JDC: Shouldn't this be larger for reaction field?
        elif explicit:
            ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
            pdb = app.PDBFile('5dfr_solv-cube_equil.pdb')
            if testName == 'pme':
                method = app.PME
                cutoff = pme_cutoff
            else:
                # Reaction field uses hard-coded cutoff
                method = app.CutoffPeriodic
                cutoff = 1.0 # nanometers; JDC: Shouldn't this be larger for reaction field?
        else:
            ff = app.ForceField('amber99sb.xml', 'amber99_obc.xml')
            pdb = app.PDBFile('5dfr_minimized.pdb')
            method = app.CutoffNonPeriodic
            cutoff = 2.0 # nanometers
        if bond_constraints == 'hbonds':
            constraints = app.HBonds
            hydrogenMass = 1.5*unit.amu
            test_parameters['constraints'] = 'HBonds'
            test_parameters['hydrogen_mass'] = '1.5'
        elif bond_constraints == 'allbonds':
            constraints = app.AllBonds # CAUTION: Constraining all bonds will perturb the effective dihedral marginal distributions.
            hydrogenMass = 4*unit.amu
            test_parameters['constraints'] = 'AllBonds'
            test_parameters['hydrogen_mass'] = '4'
        else:
            raise ValueError(f"bond_constraints must be one of 'hbonds', 'allbonds': found {bond_constraints}")
            
        test_parameters['cutoff'] = cutoff

        positions = pdb.positions
        system = ff.createSystem(pdb.topology, nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints, hydrogenMass=hydrogenMass)

    return system, positions, test_parameters

def serializeTest(directory=None, system=None, integrator=None, state=None, coredata=None, metadata=None):
    """Serialize test system for benchmarking on Folding@home.

    Parameters
    ----------
    directory : str
        Directory to write serialized test XML files
    system : openmm.System
        The System to serialize
    integrator : openmm.Integrator
        The Integrator to serialize
    state : openmm.State
        The State to serialize
    coredata : dict
        Core parameters to serialize to XML
    metadata : dict
        Metadata to dump to YAML
    """
    os.makedirs(directory, exist_ok=True)
    import openmm
    def serialize(obj, filename):
        if filename.endswith('.bz2'):
            import bz2
            with bz2.open(filename, 'wt') as outfile:
                outfile.write(openmm.XmlSerializer.serialize(obj))
        elif filename.endswith('.gz'):
            import gzip
            with gzip.open(filename, 'wt') as outfile:
                outfile.write(openmm.XmlSerializer.serialize(obj))
        else:
            with open(filename, 'wt') as outfile:
                outfile.write(openmm.XmlSerializer.serialize(obj))

    serialize(system, os.path.join(directory, 'system.xml.bz2'))
    serialize(integrator, os.path.join(directory, 'integrator.xml.bz2'))
    serialize(state, os.path.join(directory, 'state.xml.bz2'))

    with open(os.path.join(directory, 'core.xml'), 'wt') as outfile:
        outfile.write('<config>\n')
        for key, value in coredata.items():
            outfile.write(f'    <{key} v="{value}"/>\n')
        outfile.write('</config>\n')

    import yaml
    with open(os.path.join(directory, 'metadata.yaml'), 'wt') as outfile:
        outfile.write(yaml.dump(metadata))

def runOneTest(testName, options):
    """Perform a single benchmarking simulation."""

    system, positions, test_parameters = retrieveTestSystem(testName, pme_cutoff=options.pme_cutoff, bond_constraints=options.bond_constraints, polarization=options.polarization, epsilon=options.epsilon)

    # Create a copy of the basic test_parameters dict (which may be cached) to report the test results
    test_result = test_parameters.copy()
    test_result['ensemble'] = options.ensemble
    test_result['precision'] = options.precision

    explicit = (testName not in ('gbsa', 'amoebagk'))
    amoeba = (testName in ('amoebagk', 'amoebapme'))
    apoa1 = testName.startswith('apoa1')
    amber = (testName.startswith('amber'))
    
    # Create the integrator
    temperature = 300*unit.kelvin
    if explicit:
        friction = 1*(1/unit.picoseconds)
    else:
        friction = 91*(1/unit.picoseconds)
    if amoeba:
        dt = 0.002*unit.picoseconds
        if options.ensemble == 'NVE':
            integ = mm.MTSIntegrator(dt, [(0,2), (1,1)])
        else:
            integ = mm.MTSLangevinIntegrator(temperature, friction, dt, [(0,2), (1,1)])
    elif amber:
        dt = 0.004*unit.picoseconds
        if options.ensemble == 'NVE':
            integ = mm.VerletIntegrator(dt)
        else:
            integ = mm.LangevinMiddleIntegrator(temperature, friction, dt)
    else:
        if options.bond_constraints == 'hbonds':
            dt = 0.004*unit.picoseconds
            if options.ensemble == 'NVE':
                integ = mm.VerletIntegrator(dt)
            else:
                integ = mm.LangevinMiddleIntegrator(temperature, friction, dt)
        elif options.bond_constraints == 'allbonds':
            dt = 0.005*unit.picoseconds
            if options.ensemble == 'NVE':
                integ = mm.VerletIntegrator(dt)
            else:
                integ = mm.LangevinIntegrator(temperature, friction, dt)
        else:
            raise ValueError(f'Unknown bond_constraints {options.bond_constraints}')

    test_result['timestep_in_fs'] = dt.value_in_unit(unit.femtoseconds)
    properties = {}
    initialSteps = 5
    platform = mm.Platform.getPlatform(options.platform)
    if options.device is not None and 'DeviceIndex' in platform.getPropertyNames():
        properties['DeviceIndex'] = options.device
        if ',' in options.device or ' ' in options.device:
            initialSteps = 250
    if options.disable_pme_stream:
        properties['DisablePmeStream'] = 'true'
    if options.opencl_platform is not None and 'OpenCLPlatformIndex' in platform.getPropertyNames():
        properties['OpenCLPlatformIndex'] = options.opencl_platform
    if (options.precision is not None) and ('Precision' in platform.getPropertyNames()):
        properties['Precision'] = options.precision

    # Add barostat if requested
    if options.ensemble == 'NPT':
        system.addForce(mm.MonteCarloBarostat(1*unit.bar, temperature, 100))

    # Create the Context
    integ.setConstraintTolerance(1e-5)
    if len(properties) > 0:
        context = mm.Context(system, integ, platform, properties)
    else:
        context = mm.Context(system, integ, platform)

    # Store information about the Platform used by the Context
    platform = context.getPlatform()
    test_result['platform'] = platform.getName()
    test_result['platform_properties'] = { property_name : platform.getPropertyValue(context, property_name) for property_name in platform.getPropertyNames() }

    # Prepare the simulation
    context.setPositions(positions)
    if amber:
        mm.LocalEnergyMinimizer.minimize(context, 100*unit.kilojoules_per_mole/unit.nanometer)
    context.setVelocitiesToTemperature(temperature)

    if options.serialize:
        # Apply constraints on positions and velocities if needing to serialize for Folding@home
        tol = 1.0e-8
        context.applyConstraints(tol)
        context.applyVelocityConstraints(tol)
        state = context.getState(positions=True, velocities=True, energy=True, forces=True, parameters=True)

    # Time integration, ensuring we trigger kernel compilation before we start timing
    steps = 20
    while True:
        elapsed_time = timeIntegration(context, steps, initialSteps)
        if elapsed_time >= 0.5*options.seconds:
            break
        if elapsed_time < 0.5:
            steps = int(steps*1.0/elapsed_time) # Integrate enough steps to get a reasonable estimate for how many we'll need.
        else:
            steps = int(steps*options.seconds/elapsed_time)
    test_result['steps'] = steps
    test_result['elapsed_time'] = elapsed_time
    time_per_step = elapsed_time * unit.seconds / steps 
    ns_per_day = (integ.getStepSize() / time_per_step) / (unit.nanoseconds/unit.day)
    test_result['ns_per_day'] = ns_per_day

    # Serialize XML files for Folding@home benchmark if requested
    if options.serialize:
        wu_duration = 5*unit.minutes
        ncheckpoints_per_wu = 2
        nsteps_per_wu = int(wu_duration / time_per_step)

        coredata = dict()
        coredata['checkpointFreq'] = int(nsteps_per_wu / ncheckpoints_per_wu)
        coredata['numSteps'] = ncheckpoints_per_wu * coredata['checkpointFreq']
        coredata['xtcFreq'] = coredata['numSteps']
        coredata['precision'] = options.precision
        coredata['xtcAtoms'] = 'solute'
        coredata['disableCheckpointStateTests'] = 1 # disable checkpoint state tests
        coredata['DisablePmeStream'] = 0 # don't disable separate PME streams
        coredata['forceTolerance'] = 200000 # for some reason, we have to make this large (TODO: Investigate why)
        coredata['energyTolerance'] = 200000 # for some reason, we have to make this large (TODO: Investigate why)
        
        forces = { force.getName() : force for force in system.getForces() }
        NONBONDED_METHODS = { 0 : 'NoCutoff', 1 : 'CutoffNonPeriodic', 2 : 'CutoffPeriodic', 3 : 'Ewald', 4 : 'PME', 5 : 'LJPME' }
        if 'NonbondedForce' in forces:
            nonbonded_method = NONBONDED_METHODS[forces['NonbondedForce'].getNonbondedMethod()]
        else:
            nonbonded_method = 'unknown'

        metadata = {
            'name' : testName,
            'description' : f'OpenMM Benchmark Suite : {testName}',
            'precision' : options.precision,
            'num_atoms' : context.getSystem().getNumParticles(),
            'nonbonded_method' : f'{nonbonded_method}',
            'timestep' : integ.getStepSize() / unit.picoseconds,
            'integrator' : f'{integ.__class__.__name__}',
            }

        directory = os.path.join(options.serialize, f'{testName}-{options.precision}')
        serializeTest(directory=directory, system=context.getSystem(), integrator=context.getIntegrator(), state=state, coredata=coredata, metadata=metadata)

    # Clean up
    del context, integ

    printTestResult(test_result, options)
    appendTestResult(test_result=test_result, filename=options.outfile)

# Parse the command line options.

platform_speeds = { mm.Platform.getPlatform(i).getName() : mm.Platform.getPlatform(i).getSpeed() for i in range(mm.Platform.getNumPlatforms()) }
PLATFORMS = [platform for platform, speed in sorted(platform_speeds.items(), key=lambda item: item[1], reverse=True)]
TESTS = ('gbsa', 'rf', 'pme', 'apoa1rf', 'apoa1pme', 'apoa1ljpme', 'amoebagk', 'amoebapme', 'amber20-dhfr', 'amber20-cellulose', 'amber20-stmv')
ENSEMBLES = ('NVE', 'NVT', 'NPT')
BOND_CONSTRAINTS = ('hbonds', 'allbonds')
PRECISIONS = ('single', 'mixed', 'double')
POLARIZATION_MODES = ('direct', 'extrapolated', 'mutual')
STYLES = ('simple', 'table')

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Run one or more benchmarks of OpenMM",
                                epilog="""
Example: run the full suite of benchmarks for the CUDA platform, printing the results as a table

    python benchmark.py --platform=CUDA --style=table

Example: run the apoa1pme benchmark for the CPU platform with a reduced cutoff distance

    python benchmark.py --platform=CPU --test=apoa1pme --pme-cutoff=0.8

Example: run the full suite in mixed precision mode, saving the results to a YAML file

    python benchmark.py --platform=CUDA --precision=mixed --outfile=benchmark.yaml""")
parser.add_argument('--platform', dest='platform', choices=PLATFORMS, help='name of the platform to benchmark')
parser.add_argument('--test', default=','.join(TESTS), dest='test', help=f'the test to perform, or comma-separated list: {TESTS} [default: all]')
parser.add_argument('--ensemble', default='NVT', dest='ensemble', help=f'the thermodynamic ensemble to simulate: {ENSEMBLES} [default: NVT]')
parser.add_argument('--pme-cutoff', default=0.9, dest='pme_cutoff', type=float, help='direct space cutoff for PME in nm [default: 0.9]')
parser.add_argument('--seconds', default=60, dest='seconds', type=float, help='target simulation length in seconds [default: 60]')
parser.add_argument('--polarization', default='mutual', dest='polarization', choices=POLARIZATION_MODES, help=f'the polarization method for AMOEBA: {POLARIZATION_MODES} [default: mutual]')
parser.add_argument('--mutual-epsilon', default=1e-5, dest='epsilon', type=float, help='mutual induced epsilon for AMOEBA [default: 1e-5]')
parser.add_argument('--bond-constraints', default='hbonds', dest='bond_constraints', help=f'hbonds: constrain bonds to hydrogen, use 1.5*amu H mass; allbonds: constrain all bonds, use 4*amu H mass, and use larger timestep. This option is ignored for AMOEBA: {BOND_CONSTRAINTS} [default: hbonds]')
parser.add_argument('--disable-pme-stream', default=False, action='store_true', dest='disable_pme_stream', help='disable use of a separate GPU stream for PME')
parser.add_argument('--device', default=None, dest='device', help='device index for CUDA, HIP, or OpenCL')
parser.add_argument('--opencl-platform', default=None, dest='opencl_platform', help='platform index for OpenCL')
parser.add_argument('--precision', default='single', dest='precision', help=f'precision modes for CUDA, HIP, or OpenCL: {PRECISIONS} [default: single]')
parser.add_argument('--style', default='simple', dest='style', choices=STYLES, help=f'output style: {STYLES} [default: simple]')
parser.add_argument('--outfile', default=None, dest='outfile', help='output filename for benchmark logging (must end with .yaml or .json)')
parser.add_argument('--serialize', default=None, dest='serialize', help='if specified, output serialized test systems for Folding@home or other uses')
parser.add_argument('--verbose', default=False, action='store_true', dest='verbose', help='if specified, print verbose output')
args = parser.parse_args()
if args.platform is None:
    parser.error('No platform specified')

# Collect system information
system_info = dict()
import socket, platform
system_info['hostname'] = socket.gethostname()
system_info['timestamp'] = datetime.now(timezone.utc).isoformat()
system_info['openmm_version'] = mm.version.version
system_info['cpuinfo'] = cpuinfo()
system_info['cpuarch'] = platform.processor()
system_info['system'] = platform.system()
# TODO: Capture information on how many CPU threads will be used

# Attempt to get GPU info
try:
    import shutil
    import subprocess
    if shutil.which('nvidia-smi') is not None:
        cmd = 'nvidia-smi --query-gpu=driver_version,gpu_name --format=csv,noheader'
        output = subprocess.check_output(cmd, shell=True, text=True)
        system_info['nvidia_driver'], system_info['gpu'] = output.strip().split(', ')
except Exception as e:
    pass

for key, value in system_info.items():
    print(f'{key}: {value}')

if args.outfile is not None:
    # Remove output file
    if os.path.exists(args.outfile):
        os.unlink(args.outfile)
    # Write system info
    appendTestResult(args.outfile, system_info=system_info)

tests = args.test.split(',')
if not set(tests).issubset(TESTS):
    parser.error(f'Available tests: {TESTS}')

precisions = args.precision.split(',')
if args.platform == 'Reference':
    precisions = ['double']
if args.platform == 'CPU':
    precisions = ['mixed']
if not set(precisions).issubset(PRECISIONS):
    parser.error(f'Available precisions: {PRECISIONS}')

ensembles = args.ensemble.split(',')
if not set(ensembles).issubset(ENSEMBLES):
    parser.error(f'Available ensembles: {ENSEMBLES}')

bond_constraints = args.bond_constraints.split(',')
if not set(bond_constraints).issubset(BOND_CONSTRAINTS):
    parser.error(f'Available bond constraints: {BOND_CONSTRAINTS}')

# Combinatorially run all requested benchmarks, ignoring combinations that cannot be run
from openmm import OpenMMException
from itertools import product

if args.style == 'simple':
    for (test, args.bond_constraints, args.ensemble, args.precision) in product(tests, bond_constraints, ensembles, precisions):
        try:
            runOneTest(test, args)
        except OpenMMException as e:
            if args.verbose:
                print(e)
            pass
elif args.style == 'table':
    print()
    print('Test              Precision   Constraints   H mass (amu)   dt (fs)   Ensemble   Platform   ns/day')
    for (test, args.bond_constraints, args.ensemble, args.precision) in product(tests, bond_constraints, ensembles, precisions):
        try:
            runOneTest(test, args)
        except OpenMMException as e:
            if args.verbose:
                print(e)
            pass
else:
    raise ValueError(f"style {args.style} unknown; must be one of ['simple', 'table']")
