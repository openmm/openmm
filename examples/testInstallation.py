from __future__ import print_function
# First make sure OpenMM is installed.


import sys
try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
except ImportError as err:
    print("Failed to import OpenMM packages:", err.message)
    print("Make sure OpenMM is installed and the library path is set correctly.")
    sys.exit()

# Create a System for the tests.

pdb = PDBFile('input.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# List all installed platforms and compute forces with each one.

numPlatforms = Platform.getNumPlatforms()
print("There are", numPlatforms, "Platforms available:")
print()
forces = [None]*numPlatforms
for i in range(numPlatforms):
    platform = Platform.getPlatform(i)
    print(i+1, platform.getName(), end=" ")
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    try:
        simulation = Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        forces[i] = simulation.context.getState(getForces=True).getForces()
        del simulation
        print("- Successfully computed forces")
    except:
        print("- Error computing forces with", platform.getName(), "platform")

# See how well the platforms agree.

if numPlatforms > 1:
    print()
    print("Median difference in forces between platforms:")
    print()
    for i in range(numPlatforms):
        for j in range(i):
            if forces[i] is not None and forces[j] is not None:
                errors = []
                for f1, f2 in zip(forces[i], forces[j]):
                    d = f1-f2
                    error = sqrt((d[0]*d[0]+d[1]*d[1]+d[2]*d[2])/(f1[0]*f1[0]+f1[1]*f1[1]+f1[2]*f1[2]))
                    errors.append(error)
                print("{} vs. {}: {:g}".format(Platform.getPlatform(j).getName(),
                                              Platform.getPlatform(i).getName(),
                                              sorted(errors)[len(errors)//2]))
