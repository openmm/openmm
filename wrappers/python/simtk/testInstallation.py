from __future__ import print_function
import os
import sys
# First make sure OpenMM is installed.

class TestingError(Exception):
    """
    An error is encountered when 
    """

def run_tests():
    try:
        from simtk.openmm.app import *
        from simtk.openmm import *
        from simtk.unit import *
    except ImportError as err:
        raise TestingError('Failed to import OpenMM packages; Make sure OpenMM\n'
                           'is installed and the library path is set correctly.\n'
                           'Error message: %s' % err.message)
    
    # Create a System for the tests.
    data_dir = os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'openmm', 'app', 'data')
    pdb = PDBFile(os.path.join(data_dir, 'test.pdb'))
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    
    # List all installed platforms and compute forces with each one.
    
    numPlatforms = Platform.getNumPlatforms()
    print("There are", numPlatforms, "Platforms available:")
    print()
    forces = [None]*numPlatforms
    platformErrors = {}
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
            platformErrors[platform.getName()] = sys.exc_info()[1]
    
    # Give details of any errors.
    
    for platform in platformErrors:
        print()
        print("%s platform error: %s" % (platform, platformErrors[platform]))
    
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

if __name__ == '__main__':

    try:
        run_tests()
    except TestingError as err:
        print('An error was encountered during testing:')
        print(err)
        sys.exit(1)
    except BaseException as err:
        print('An unexpected error was encountered during testing:')
        print(err)
        sys.exit(1)
