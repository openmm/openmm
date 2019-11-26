import math
import sys

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

def VerletIntegrator2(timeStep):

    integrator = CustomIntegrator(timeStep)

    integrator.addUpdateContextState()
    integrator.addComputePerDof('v', 'v+0.5*dt*f/m')
    integrator.addComputePerDof('x', 'x+dt*v')
    integrator.addComputePerDof('v', 'v+0.5*dt*f/m')

    return integrator

def MultiStepVerletIntegrator2(timeStep, slowPeriod):

    integrator = CustomIntegrator(timeStep)

    integrator.addGlobalVariable('slowStep', 0)
    integrator.addPerDofVariable('fe', 0)

    integrator.addComputeGlobal('slowStep', 'slowStep + 1')
    integrator.beginIfBlock(f'slowStep = {slowPeriod}')
    integrator.addComputeGlobal('slowStep', '0')
    integrator.endBlock()

    integrator.addUpdateContextState()

    integrator.addComputePerDof('fe', 'f0')
    integrator.beginIfBlock('slowStep = 0')
    integrator.addComputePerDof('fe', 'fe + f1')
    integrator.endBlock()

    integrator.addComputePerDof('v', 'v+0.5*dt*fe/m');
    integrator.addComputePerDof('x', 'x+dt*v');

    integrator.addComputePerDof('fe', 'f0')
    integrator.beginIfBlock('slowStep = 0')
    integrator.addComputePerDof('fe', 'fe + f1')
    integrator.endBlock()

    integrator.addComputePerDof('v', 'v+0.5*dt*fe/m');

    return integrator

def setForceGroups(system):

    for force in system.getForces():
        force.setForceGroup(0) # 0x00000001 = 1
        if isinstance(force, NonbondedForce):
            force.setReciprocalSpaceForceGroup(1) # 0x00000010 = 2


if __name__ == '__main__':

    integratorType = sys.argv[1]

    prmtop = AmberPrmtopFile('input.prmtop')
    system = prmtop.createSystem(nonbondedMethod=PME, rigidWater=False)

    if integratorType == 'v1':
        integrator = VerletIntegrator(0.5*femtoseconds)
    elif integratorType == 'v2':
        integrator = VerletIntegrator2(0.5*femtoseconds)
    elif integratorType == 'l1':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator(0.5*femtoseconds, 1, 2, 1)
    elif integratorType == 'l2':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator2(0.5*femtoseconds, 1)
    elif integratorType == 'm1':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator(0.5*femtoseconds, 1, 2, 2)
    elif integratorType == 'm2':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator2(0.5*femtoseconds, 2)

    platform = Platform.getPlatformByName('CUDA')
    #properties = {'DeterministicForces': 'true', 'Precision': 'double'}
    properties = {}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    inpcrd = AmberInpcrdFile('input.inpcrd')
    simulation.context.setPositions(inpcrd.positions)

    simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, speed=True))
    simulation.step(10000)