from __future__ import print_function
import openmm.app as app
import openmm as mm
import openmm.unit as unit
from datetime import datetime
import os
from argparse import ArgumentParser

def cpuinfo():
    """Return CPU info"""
    import os, platform, subprocess, re
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

    import os
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
    elif options.style == 'rich':
        options.rich_table.add_row(
            f"[red]{test_result['test']}",
            f"[blue]{test_result['precision']}",
            f"[blue_violet]{test_result['constraints']}",
            f"[blue_violet]{test_result['hydrogen_mass']}",
            f"[blue_violet]{test_result['timestep_in_fs']:.1f}",
            f"[orange1]{test_result['ensemble']}",
            f"[blue]{test_result['platform']}",
            f"[bold][green]{test_result['ns_per_day']:.1f}",
        )
        options.rich_live.refresh()
    else:
        raise ValueError(f"style '{style}' must be one of ['legacy', 'rich']")

def timeIntegration(context, steps, initialSteps):
    """Integrate a Context for a specified number of steps, then return how many seconds it took."""
    context.getIntegrator().step(initialSteps) # Make sure everything is fully initialized
    context.getState(getEnergy=True)
    start = datetime.now()
    context.getIntegrator().step(steps)
    context.getState(getEnergy=True)
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
@functools.cache
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

    """
    explicit = (testName not in ('gbsa', 'amoebagk'))
    amoeba = (testName in ('amoebagk', 'amoebapme'))
    apoa1 = testName.startswith('apoa1')
    amber = (testName.startswith('amber'))
    hydrogenMass = None

    # Create the System.
    if amoeba:
        constraints = None
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
    elif amber:
        dirname = downloadAmberSuite()
        names = {'amber20-dhfr':'JAC', 'amber20-factorix':'FactorIX', 'amber20-cellulose':'Cellulose', 'amber20-stmv':'STMV'}
        fileName = names[testName]
        prmtop = app.AmberPrmtopFile(os.path.join(dirname, f'PME/Topologies/{fileName}.prmtop'))
        inpcrd = app.AmberInpcrdFile(os.path.join(dirname, f'PME/Coordinates/{fileName}.inpcrd'))
        positions = inpcrd.positions
        method = app.PME
        constraints = app.HBonds
        hydrogenMass = 1.5*unit.amu
        cutoff = pme_cutoff
        system = prmtop.createSystem(nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints)
        if inpcrd.boxVectors is not None:
            system.setDefaultPeriodicBoxVectors(*inpcrd.boxVectors)
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
                cutoff = 1*unit.nanometers # JDC: Shouldn't this be larger for reaction field?
        elif explicit:
            ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
            pdb = app.PDBFile('5dfr_solv-cube_equil.pdb')
            if testName == 'pme':
                method = app.PME
                cutoff = pme_cutoff
            else:
                # Reaction field uses hard-coded cutoff
                method = app.CutoffPeriodic
                cutoff = 1*unit.nanometers # JDC: Shouldn't this be larger for reaction field?
        else:
            ff = app.ForceField('amber99sb.xml', 'amber99_obc.xml')
            pdb = app.PDBFile('5dfr_minimized.pdb')
            method = app.CutoffNonPeriodic # JDC: Shouldn't this be app.NoCutoff?
            cutoff = 2*unit.nanometers
        if bond_constraints == 'hbonds':
            constraints = app.HBonds
            hydrogenMass = 1.5*unit.amu
        elif bond_constraints == 'allbonds':
            constraints = app.AllBonds # JDC: Isn't constraining all bonds problematic?
            hydrogenMass = 4*unit.amu
        else:
            raise ValueError(f"bond_constraints must be one of 'hbonds', 'allbonds': found {bond_constraints}")

        positions = pdb.positions
        system = ff.createSystem(pdb.topology, nonbondedMethod=method, nonbondedCutoff=cutoff, constraints=constraints, hydrogenMass=hydrogenMass)

    return system, positions

def runOneTest(testName, options):
    """Perform a single benchmarking simulation."""

    # Return if we do not support the requested precision
    if ((options.platform == 'Reference') and (options.precision != 'double')) \
       or ((options.platform == 'CPU') and (options.precision != 'mixed')):
        return

    system, positions = retrieveTestSystem(testName, pme_cutoff=options.pme_cutoff, bond_constraints=options.bond_constraints, polarization=options.polarization, epsilon=options.epsilon)

    test_result = dict()
    test_result['test'] = testName

    explicit = (testName not in ('gbsa', 'amoebagk'))
    amoeba = (testName in ('amoebagk', 'amoebapme'))
    apoa1 = testName.startswith('apoa1')
    amber = (testName.startswith('amber'))
    hydrogenMass = None

    if amoeba:
        test_result['epsilon'] = options.epsilon
    elif (testName in ['pme', 'ljpme']) or testName.startswith('amber'):
        test_result['cutoff'] = options.pme_cutoff

    if amoeba:
        test_result['constraints'] = 'None'
        test_result['hydrogen_mass'] = '1'
    elif options.bond_constraints == 'hbonds':
        test_result['constraints'] = 'HBonds'
        test_result['hydrogen_mass'] = '1.5'
    elif options.bond_constraints == 'allbonds':
        test_result['constraints'] = 'AllBonds'
        test_result['hydrogen_mass'] = '4'
    else:
        raise ValueError(f'Unknown bond_constraints: {bond_constraints}')

    test_result['ensemble'] = options.ensemble
    platform = mm.Platform.getPlatformByName(options.platform)
    test_result['platform'] = options.platform
    test_result['precision'] = options.precision

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
            hydrogenMass = 4*unit.amu
            if options.ensemble == 'NVE':
                integ = mm.VerletIntegrator(dt)
            else:
                integ = mm.LangevinIntegrator(temperature, friction, dt)
        else:
            raise ValueError(f'Unknown bond_constraints {options.bond_constraints}')

    test_result['timestep_in_fs'] = dt.value_in_unit(unit.femtoseconds)
    properties = {}
    initialSteps = 5
    if options.device is not None and platform.getName() in ('CUDA', 'OpenCL'):
        properties['DeviceIndex'] = options.device
        if ',' in options.device or ' ' in options.device:
            initialSteps = 250
    if options.precision is not None and platform.getName() in ('CUDA', 'OpenCL'):
        properties['Precision'] = options.precision

    # Add barostat if requested
    if options.ensemble == 'NPT':
        system.addForce(mm.MonteCarloBarostat(1*unit.bar, temperature, 100))

    # Run the simulation.

    integ.setConstraintTolerance(1e-5)
    if len(properties) > 0:
        context = mm.Context(system, integ, platform, properties)
    else:
        context = mm.Context(system, integ, platform)
    context.setPositions(positions)
    if amber:
        mm.LocalEnergyMinimizer.minimize(context, 100*unit.kilojoules_per_mole/unit.nanometer)
    context.setVelocitiesToTemperature(temperature)

    steps = 20
    while True:
        time = timeIntegration(context, steps, initialSteps)
        if time >= 0.5*options.seconds:
            break
        if time < 0.5:
            steps = int(steps*1.0/time) # Integrate enough steps to get a reasonable estimate for how many we'll need.
        else:
            steps = int(steps*options.seconds/time)
    test_result['steps'] = steps
    test_result['elapsed_time'] = time
    ns_per_day = (dt*steps*86400/time).value_in_unit(unit.nanoseconds)
    test_result['ns_per_day'] = ns_per_day

    # Clean up
    del context, integ

    printTestResult(test_result, options)
    appendTestResult(test_result=test_result, filename=options.outfile)

# Parse the command line options.

platform_speeds = { mm.Platform.getPlatform(i).getName() : mm.Platform.getPlatform(i).getSpeed() for i in range(mm.Platform.getNumPlatforms()) }
PLATFORMS = [platform for platform, speed in sorted(platform_speeds.items(), key=lambda item: item[1], reverse=True)]
TESTS = ('gbsa', 'rf', 'pme', 'apoa1rf', 'apoa1pme', 'apoa1ljpme', 'amoebagk', 'amoebapme', 'amber20-dhfr',  'amber20-factorix', 'amber20-cellulose', 'amber20-stmv')
ENSEMBLES = ('NVE', 'NVT', 'NPT')
BOND_CONSTRAINTS = ('hbonds', 'allbonds')
PRECISIONS = ('single', 'mixed', 'double')
POLARIZATION_MODES = ('direct', 'extrapolated', 'mutual')

parser = ArgumentParser()
parser.add_argument('--platform', default=None, dest='platform', help=f'name of the platform to benchmark, or comma-separated list (CUDA,OpenCL): {PLATFORMS} [default: all]')
parser.add_argument('--test', default=None, dest='test', help=f'the test to perform, or comma-separated list (pme,rf): {TESTS} [default: all]')
parser.add_argument('--ensemble', default=None, dest='ensemble', help=f'the thermodynamic ensemble to simulate: {ENSEMBLES} [default: all]')
parser.add_argument('--pme-cutoff', default=0.9, dest='pme_cutoff', type=float, help='direct space cutoff for PME in nm [default: 0.9]')
parser.add_argument('--seconds', default=60, dest='seconds', type=float, help='target simulation length in seconds [default: 60]')
parser.add_argument('--polarization', default='mutual', dest='polarization', choices=POLARIZATION_MODES, help='the polarization method for AMOEBA: {POLARIZATION_MODES} [default: mutual]')
parser.add_argument('--mutual-epsilon', default=1e-5, dest='epsilon', type=float, help='mutual induced epsilon for AMOEBA [default: 1e-5]')
parser.add_argument('--bond-constraints', default=None, dest='bond_constraints', help='hbonds: constrain bonds to hydrogen, use 1.5*amu H mass; allbonds: constrain all bonds, use 4*amu H mass [default: all]')
parser.add_argument('--device', default=None, dest='device', help='device index for CUDA or OpenCL')
parser.add_argument('--precision', default=None, dest='precision', help=f'precision mode for CUDA or OpenCL: {PRECISIONS} [default: all]')
parser.add_argument('--style', default='simple', dest='style', help='output style [default: simple]')
parser.add_argument('--outfile', default=None, dest='outfile', help='output filename for benchmark logging (must end with .yaml or .json)')
args = parser.parse_args()

# Collect system information
system_info = dict()
import socket, platform
system_info['hostname'] = socket.gethostname()
system_info['timestamp'] = datetime.utcnow().isoformat()
system_info['openmm_version'] = mm.version.version
system_info['cpuinfo'] = cpuinfo()
system_info['cpuarch'] = platform.processor()
system_info['system'] = platform.system()
# TODO: Capture information on how many CPU threads will be used

# Attempt to get GPU info
try:
    import subprocess
    cmd = 'nvidia-smi --query-gpu=driver_version,gpu_name --format=csv,noheader'
    output = subprocess.check_output(cmd, shell=True, text=True)
    system_info['nvidia_driver'], system_info['gpu'] = output.strip().split(', ')
except Exception as e:
    pass

for key, value in system_info.items():
    print(f'{key}: {value}')

if args.outfile is not None:
    # Remove output file
    import os
    if os.path.exists(args.outfile):
        os.unlink(args.outfile)
    # Write system info
    appendTestResult(args.outfile, system_info=system_info)

if args.platform is not None:
    # Use specified subset of platforms
    requested_platforms = args.platform.split(',')
    if not set(requested_platforms).issubset(PLATFORMS):
        parser.error(f'Available platforms: {PLATFORMS}')
    PLATFORMS = requested_platforms

if args.test is not None:
    # Use specified subset of tests
    requested_tests = args.test.split(',')
    if not set(requested_tests).issubset(TESTS):
        parser.error(f'Available tests: {TESTS}')
    TESTS = requested_tests

if args.precision is not None:
    # Use specified subset of precisions
    requested_precisions = args.precision.split(',')
    if not set(requested_precisions).issubset(PRECISIONS):
        parser.error(f'Available precisions: {PRECISIONS}')
    PRECISIONS = requested_precisions

if args.ensemble is not None:
    # Use specified subset of ensembles
    requested_ensembles = args.ensemble.split(',')
    if not set(requested_ensembles).issubset(ENSEMBLES):
        parser.error(f'Available ensembles: {ENSEMBLES}')
    ENSEMBLES = requested_ensembles

if args.bond_constraints is not None:
    # Use specified subset of constraints
    requested_bond_constraints = args.bond_constraints.split(',')
    if not set(requested_bond_constraints).issubset(BOND_CONSTRAINTS):
        parser.error(f'Available bond constraints: {BOND_CONSTRAINTS}')
    BOND_CONSTRAINTS = requested_bond_constraints

# Combinatorially run all requested benchmarks, ignoring combinations that cannot be run
from openmm import OpenMMException
from itertools import product

if args.style == 'simple':
    for (test, args.bond_constraints, args.ensemble, args.platform, args.precision) in product(TESTS, BOND_CONSTRAINTS, ENSEMBLES, PLATFORMS, PRECISIONS):
        try:
            runOneTest(test, args)
        except OpenMMException as e:
            pass
elif args.style == 'rich':
    from rich.live import Live
    from rich.table import Table

    from rich import box
    table = Table(box=box.SIMPLE)
    table.add_column("[red]Test", width=24)
    table.add_column("[blue]Precision")
    table.add_column("[blue_violet]Constraints")
    table.add_column("[blue_violet]H mass (amu)")
    table.add_column("[blue_violet]dt (fs)")
    table.add_column("[blue]Platform", width=9)
    table.add_column("[orange1]Ensemble")
    table.add_column("[green]ns/day", justify='right', width=10)
    setattr(args, 'rich_table', table)

    with Live(table, auto_refresh=False, vertical_overflow='visible') as live:
        for (test, args.bond_constraints, args.ensemble, args.platform, args.precision) in product(TESTS, BOND_CONSTRAINTS, ENSEMBLES, PLATFORMS, PRECISIONS):
            try:
                setattr(args, 'rich_live', live)
                runOneTest(test, args)
            except OpenMMException as e:
                pass
else:
    raise ValueError(f"style {args.style} unkown; must be one of ['simple', 'rich']")
