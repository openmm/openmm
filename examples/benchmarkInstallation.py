from __future__ import print_function
import sys
import time
import bz2

import sys
try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
except ImportError as err:
    print("Failed to import OpenMM packages:", err.message)
    print("Make sure OpenMM is installed and the library path is set correctly.")
    sys.exit()

N_STEPS = 5000
nonbondedMethod = PME
nonbondedCutoff = 1*nanometer
constraints = HBonds
timestep = 2*femtosecond
temperature = 300*kelvin

pdb = PDBFile(bz2.BZ2File('dhfr.pdb.bz2'))
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod,
                                 nonbondedCutoff=nonbondedCutoff, constraints=constraints)

print('System Setup')
print('  NonbondedMethod: %s' % nonbondedMethod)
print('  NonbondedCutoff: %s' % nonbondedCutoff)
print('  Constraints:     %s' % constraints)
print('  Timestep:        %s' % timestep)
print('  Number of atoms: %d' % pdb.topology._numAtoms)
print('')

for i in range(Platform.getNumPlatforms())[::-1]:
    platform = Platform.getPlatform(i)
    integrator = LangevinIntegrator(temperature, 1/picosecond, timestep)
    context = Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    context.setVelocitiesToTemperature(temperature)

    print('Platform: %s' % platform.getName())
    for key in platform.getPropertyNames():
        print('  %s: %s' % (key, platform.getPropertyValue(context, key)))

    startTime = time.time()
    integrator.step(N_STEPS)
    endTime = time.time()
    wallTime = ((endTime - startTime)*seconds).value_in_unit(days)
    simTime = context.getState().getTime().value_in_unit(nanoseconds)
    
    print()
    print('  Performance = %.3f ns/day' % (simTime / wallTime))
    print()

