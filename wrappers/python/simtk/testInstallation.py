from __future__ import print_function
from __future__ import absolute_import
import os
import sys
# First make sure OpenMM is installed.

try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
except ImportError as err:
    simtk_import_failed = True
    simtk_import_error = str(err)
else:
    simtk_import_failed = False


def run_tests():
    """
    Runs a set of tests to determine which platforms are available and tests the
    relative accuracy between them. This can be used to determine if the Python
    API is installed and working properly, as well as the fidelity of the
    underlying OpenMM libraries with respect to computing energies and forces on
    the different platforms supported by your installation.

    This test prints the available platforms and the relative force errors
    between them for a test system. If a problem is detected, TestingError is
    raised.
    """
    print()
    print('OpenMM Version:', Platform.getOpenMMVersion())
    print('Git Revision:', version.git_revision)
    print()

    # Create a System for the tests.
    data_dir = os.path.join(os.path.abspath(os.path.split(__file__)[0]), 'openmm', 'app', 'data')
    pdb = PDBFile(os.path.join(data_dir, 'test.pdb'))
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=LJPME, nonbondedCutoff=1*nanometer, constraints=HBonds, ewaldErrorTolerance=1e-4)

    # List all installed platforms and compute forces with each one.

    numPlatforms = Platform.getNumPlatforms()
    print("There are", numPlatforms, "Platforms available:")
    print()
    forces = [None] * numPlatforms
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

    errorsOk = True
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
                    median = sorted(errors)[len(errors)//2]
                    if median < 1e-4:
                        message = ''
                    else:
                        message = '  *** LARGE DIFFERENCE **'
                        errorsOk = False
                    print("{0} vs. {1}: {2:g}{3}".format(Platform.getPlatform(j).getName(),
                                                  Platform.getPlatform(i).getName(), median, message))
    if errorsOk:
        print()
        print('All differences are within tolerance.')


def main():
    if simtk_import_failed:
        print('Failed to import OpenMM packages; OpenMM will not work.\n'            
              'Make sure OpenMM is installed and the library path is set correctly.'
              '\n\nError message: %s' % simtk_import_error,
              file=sys.stderr)
        sys.exit(1)

    try:
        run_tests()
    except Exception as err:
        print('Problem with OpenMM installation '
              'encountered. OpenMM will not work until the problem '
              'has been fixed.\n\n',
              file=sys.stderr)
        print('Error message: %s' % str(err), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
