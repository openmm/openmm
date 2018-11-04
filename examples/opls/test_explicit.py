from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import stdout
import os, time, shutil
from datetime import datetime

print("Started at: " + str(time.asctime()))
start=datetime.now()

shutil.copyfile('als_peptide_opls_rst1.dms','als_peptide_opls_rst2.dms')
path = os.path.join(os.path.dirname(__file__), 'als_peptide_opls_rst2.dms')
testDes = DesmondDMSFile(path)
system = testDes.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*nanometer, OPLS=True)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('OpenCL')
simulation = Simulation(testDes.topology, system, integrator,platform)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)
print "Using platform %s" % simulation.context.getPlatform().getName()
state = simulation.context.getState(getEnergy = True)
print(state.getPotentialEnergy())
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(DCDReporter('als_peptide_opls.dcd', 100))
#simulation.minimizeEnergy(maxIterations=1000)
simulation.step(10000)
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
#print "Updating positions and velocities"
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
