import copy
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

def MultiStepVerletIntegrator4(timeStep):

    integrator = CustomIntegrator(timeStep)

    integrator.addPerDofVariable('fe', 0)

    integrator.addUpdateContextState()

    integrator.addComputePerDof('v', 'v+0.25*dt*f0/m')
    integrator.addComputePerDof('x', 'x+0.5*dt*v')
    integrator.addComputePerDof('v', 'v+0.5*dt*f0/m')

    integrator.addUpdateContextState()

    # integrator.addComputePerDof('fe', 'f0 + 2*f1') # This doesn't work
    integrator.addComputePerDof('fe', 'f0')
    integrator.addComputePerDof('fe', 'fe + f1') # Does the same as MultiStepVerletIntegrator3
    # integrator.addComputePerDof('fe', 'fe + 2*f1') # This is what it has to do

    integrator.addComputePerDof('v', 'v+0.25*dt*fe/m')
    integrator.addComputePerDof('x', 'x+0.5*dt*v')

    # integrator.addComputePerDof('fe', 'f0 + 2*f1') # This doesn't work
    integrator.addComputePerDof('fe', 'f0')
    integrator.addComputePerDof('fe', 'fe + f1') # Does the same as MultiStepVerletIntegrator3
    # integrator.addComputePerDof('fe', 'fe + 2*f1') # This is what it has to do

    integrator.addComputePerDof('v', 'v+0.25*dt*fe/m')

    return integrator

def setForceGroups(system):

    for force in system.getForces():
        force.setForceGroup(0) # 0b00000001 = 1
        if isinstance(force, NonbondedForce):
            force.setReciprocalSpaceForceGroup(1) # 0b00000010 = 2

def addScaledForce(system):

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            new_force = copy.deepcopy(force)

            # Scale PME by 2, i.e. scale charges by sqrt(2)
            for i in range(new_force.getNumParticles()):
                charge, sigma, epsilon = new_force.getParticleParameters(i)
                charge *= math.sqrt(2)
                new_force.setParticleParameters(i, charge, sigma, epsilon)

            new_force.setForceGroup(2) # 0b00000100 = 4 <-- not used
            new_force.setReciprocalSpaceForceGroup(3) # 0b00001000 = 8 <-- slowForce

            system.addForce(new_force)

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
    elif integratorType == 'm3':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator3(1.0*femtoseconds, 1, 2)
    elif integratorType == 'm4':
        setForceGroups(system)
        integrator = MultiStepVerletIntegrator4(1.0*femtoseconds)
    elif integratorType == 'm5':
        setForceGroups(system)
        addScaledForce(system)
        integrator = MultiStepVerletIntegrator3(1.0*femtoseconds, 1, 8)
    elif integratorType == 'm6':
        setForceGroups(system)
        addScaledForce(system)
        integrator = MultiStepVerletIntegrator3(1.0*femtoseconds, 1, 8)

    platform = Platform.getPlatformByName('CUDA')
    #properties = {'DeterministicForces': 'true', 'Precision': 'double'}
    properties = {}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

    if integratorType == 'm6':
        simulation.context.setDefaulfForceGroups(3) # 0b00000011 = 3 <-- fastForces + slowForces

    inpcrd = AmberInpcrdFile('input.inpcrd')
    simulation.context.setPositions(inpcrd.positions)

    simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, speed=True))
    simulation.step(10000)