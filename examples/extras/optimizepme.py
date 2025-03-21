"""
optimizepme.py: Optimizes parameters for PME simulations

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970.

Portions copyright (c) 2013-2025 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import itertools
import math
from datetime import datetime

def optimizePME(system, integrator, positions, platform, properties, minCutoff, maxCutoff):
    """Run a series of simulations using different parameters to see which give the best performance.
    
    When running a simulation with PME, different combinations of parameters may give equivalent accuracy
    but differ in performance.  In particular:
    
    1. The nonbonded cutoff does not affect the accuracy of the Coulomb interaction with PME.  You can
    freely vary the cutoff distance, and OpenMM will automatically select internal parameters to give
    whatever accuracy has been selected with the ewaldErrorTolerance parameter.  (The cutoff does affect
    other nonbonded interactions, such as Lennard-Jones, so this generally places a lower limit on the
    cutoffs you consider acceptable.)
    2. In some cases, OpenMM can perform reciprocal space calculations on the CPU at the same time it is
    doing direct space calculations on the GPU.  Depending on your hardware, this might or might not
    be faster.
    
    This function runs a series of simulations to measure the performance of simulating a particular system
    on the current hardware.  This allows you to choose the combination of parameters that give the
    best performance while still providing the required accuracy.  The function prints out the results of
    each simulation, along with a final recommendation of the best parameters to use.  On exit, the
    system and properties arguments will have been modified to use the recommended parameters.
    
    Parameters:
     - system (System) the System to simulate
     - integrator (Integrator) the Integrator to use for simulating it
     - positions (list) the initial particle positions
     - platform (Platform) the Platform to use for running the simulation
     - properties (dict) any platform-specific properties you want to specify
     - minCutoff (distance) the minimum cutoff distance to try
     - maxCutoff (distance) the maximum cutoff distance to try
    """
    if unit.is_quantity(minCutoff):
        minCutoff = minCutoff.value_in_unit(unit.nanometers)
    if unit.is_quantity(maxCutoff):
        maxCutoff = maxCutoff.value_in_unit(unit.nanometers)
    
    # Find the NonbondedForce or AmoebaMultipoleForce to optimize.
    
    nonbonded = None
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            nonbonded = force
            if nonbonded.getNonbondedMethod() != mm.NonbondedForce.PME:
                raise ValueError('The System does not use PME')
            break
        if isinstance(force, mm.AmoebaMultipoleForce):
            nonbonded = force
            if nonbonded.getNonbondedMethod() != mm.AmoebaMultipoleForce.PME:
                raise ValueError('The System does not use PME')
            nonbonded.setAEwald(0)
            break
    if nonbonded is None:
        raise ValueError('The System does not include a NonbondedForce or AmoebaMultipoleForce')
    errorTolerance = nonbonded.getEwaldErrorTolerance()
    canUseCpuPme = (isinstance(nonbonded, mm.NonbondedForce) and platform.supportsKernels(['CalcPmeReciprocalForce']))
    if platform.getName() == 'CUDA':
        cpuPmeProperty = 'CudaUseCpuPme'
    else:
        cpuPmeProperty = 'OpenCLUseCpuPme'

    # Build a list of cutoff distances to try.
    
    gpuCutoffs = set()
    cpuCutoffs = set()
    gpuCutoffs.add(minCutoff)
    cpuCutoffs.add(minCutoff)
    vec1, vec2, vec3 = system.getDefaultPeriodicBoxVectors()
    errorTolerance5 = math.pow(errorTolerance, 0.2)
    boxDimensions = [x.value_in_unit(unit.nanometers) for x in (vec1[0], vec2[1], vec3[2])]
    for boxSize in boxDimensions: # Loop over the three dimensions of the periodic box.
        for gridSize in itertools.count(start=5, step=1): # Loop over possible sizes of the PME grid.
            # Determine whether this is a legal size for the FFT.
            
            unfactored = gridSize
            for factor in (2, 3, 5, 7):
                while unfactored > 1 and unfactored%factor == 0:
                    unfactored /= factor
            if unfactored not in (1, 11, 13):
                continue
            
            # Compute the smallest cutoff that will give this grid size.
            
            alpha = 1.5*gridSize*errorTolerance5/boxSize
            cutoff = math.sqrt(-math.log(2*errorTolerance))/alpha
            cutoff = 0.001*int(cutoff*1000) # Round up to the next picometer to avoid roundoff errors.
            if cutoff < minCutoff:
                break
            if cutoff < maxCutoff:
                cpuCutoffs.add(cutoff)
                if unfactored == 1:
                    gpuCutoffs.add(cutoff)
    gpuCutoffs = sorted(gpuCutoffs)
    cpuCutoffs = sorted(cpuCutoffs)
    
    # Select a length for the simulations so they will each take about 10 seconds.
    
    print()
    print('Selecting a length for the test simulations... ')
    nonbonded.setCutoffDistance(math.sqrt(minCutoff*maxCutoff))
    properties[cpuPmeProperty] = 'false'
    context = _createContext(system, integrator, positions, platform, properties)
    steps = 20
    time = 0.0
    while time < 8.0 or time > 12.0:
        time = _timeIntegrator(context, steps)
        steps = int(steps*10.0/time)
    print(steps, 'steps')
    
    # Run the simulations.
    
    print()
    print('Running simulations with standard PME')
    print()
    results = []
    properties[cpuPmeProperty] = 'false'
    gpuTimes = _timeWithCutoffs(system, integrator, positions, platform, properties, nonbonded, gpuCutoffs, steps)
    for time, cutoff in zip(gpuTimes, gpuCutoffs):
        results.append((time, cutoff, 'false'))
    if canUseCpuPme:
        print()
        print('Running simulations with CPU based PME')
        print()
        properties[cpuPmeProperty] = 'true'
        cpuTimes = _timeWithCutoffs(system, integrator, positions, platform, properties, nonbonded, cpuCutoffs, steps)
        for time, cutoff in zip(cpuTimes, cpuCutoffs):
            results.append((time, cutoff, 'true'))
    
    # Rerun the fastest configurations to make sure the results are consistent.
    
    print()
    print('Confirming results for best configurations')
    print()
    results.sort(key=lambda x: x[0])
    finalResults = []
    for time, cutoff, useCpu in results[:5]:
        nonbonded.setCutoffDistance(cutoff)
        properties[cpuPmeProperty] = useCpu
        context = _createContext(system, integrator, positions, platform, properties)
        time2 = _timeIntegrator(context, steps)
        time3 = _timeIntegrator(context, steps)
        medianTime = sorted((time, time2, time3))[1]
        finalResults.append((medianTime, cutoff, useCpu))
        print('Cutoff=%g, %s=%s' % (cutoff, cpuPmeProperty, useCpu))
        print('Times: %g, %g, %g' % (time, time2, time3))
        print('Median time: %g' % medianTime)
        print()
    
    # Select the best configuration.
    
    finalResults.sort(key=lambda x: x[0])
    best = finalResults[0]
    nonbonded.setCutoffDistance(best[1])
    properties[cpuPmeProperty] = best[2]
    print('Best configuration:')
    print()
    print('Cutoff=%g nm, %s=%s' % (best[1], cpuPmeProperty, best[2]))
    print()


def _createContext(system, integrator, positions, platform, properties):
    integrator = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(integrator))
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    return context


def _timeIntegrator(context, steps):
    context.getIntegrator().step(5) # Make sure everything is fully initialized
    context.getState(getEnergy=True)
    start = datetime.now()
    context.getIntegrator().step(steps)
    context.getState(getEnergy=True)
    end = datetime.now()
    return (end-start).total_seconds()


def _timeWithCutoffs(system, integrator, positions, platform, properties, nonbonded, cutoffs, steps):
    times = []
    for cutoff in cutoffs:
        nonbonded.setCutoffDistance(cutoff)
        context = _createContext(system, integrator, positions, platform, properties)
        time = _timeIntegrator(context, steps)
        print('cutoff=%g, time=%g' % (cutoff, time))
        times.append(time)
        if len(times) > 3 and times[-1] > times[-2] > times[-3] > times[-4]:
            # It's steadily getting slower as we increase the cutoff, so stop now.
            break
    return times
