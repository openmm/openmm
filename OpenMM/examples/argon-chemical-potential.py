# Example illustrating use of Free Energy OpenMM plugin.
# Calculation of chemical potential of argon in a periodic box.
#
# AUTHOR
#
# John D. Chodera <jchodera@berkeley.edu>
#
# REQUIREMENTS
#
# numpy - Scientific computing package - http://numpy.scipy.org
# pymbar - Python implementation of MBAR - https://simtk.org/home/pymbar
#
# REFERENCES
#
# [1] Michael R. Shirts and John D. Chodera. Statistically optimal analysis of samples from multiple equilibrium states.
# J. Chem. Phys. 129:124105 (2008)  http://dx.doi.org/10.1063/1.2978177

from __future__ import print_function
from simtk.openmm import *
from simtk.unit import *
import numpy

# =============================================================================
# Specify simulation parameters
# =============================================================================

nparticles = 216 # number of particles

mass = 39.9 * amu # mass
sigma = 3.4 * angstrom # Lennard-Jones sigma
epsilon = 0.238 * kilocalories_per_mole # Lennard-Jones well-depth
charge = 0.0 * elementary_charge # argon model has no charge

nequil_steps = 5000 # number of dynamics steps for equilibration
nprod_steps = 500 # number of dynamics steps per production iteration
nprod_iterations = 50 # number of production iterations per lambda value

reduced_density = 0.85 # reduced density rho*sigma^3
temperature = 300 * kelvin # temperature
pressure = 1 * atmosphere # pressure
collision_rate = 5.0 / picosecond # collision rate for Langevin thermostat
timestep = 1 * femtosecond # integrator timestep

# Specify lambda values for alchemically-modified intermediates.
lambda_values = [0.0, 0.25, 0.5, 0.75, 1.0]
nlambda = len(lambda_values)

# =============================================================================
# Compute box size.
# =============================================================================

volume = nparticles*(sigma**3)/reduced_density
box_edge = volume**(1.0/3.0)
cutoff = min(box_edge*0.49, 2.5*sigma) # Compute cutoff
print("sigma = %s" % sigma)
print("box_edge = %s" % box_edge)
print("cutoff = %s" % cutoff)

# =============================================================================
# Build systems at each alchemical lambda value.
# =============================================================================

print("Building alchemically-modified systems...")
alchemical_systems = list() # alchemical_systems[i] is the alchemically-modified System object corresponding to lambda_values[i]
for lambda_value in lambda_values:
   # Create argon system where first particle is alchemically modified by lambda_value.
   system = System()
   system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))
   force = CustomNonbondedForce("""4*epsilon*l12*(1/((alphaLJ*(1-l12) + (r/sigma)^6)^2) - 1/(alphaLJ*(1-l12) + (r/sigma)^6));
                                   sigma=0.5*(sigma1+sigma2);
                                   epsilon=sqrt(epsilon1*epsilon2);
                                   alphaLJ=0.5;
                                   l12=1-(1-lambda)*step(useLambda1+useLambda2-0.5)""")
   force.addPerParticleParameter("sigma")
   force.addPerParticleParameter("epsilon")
   force.addPerParticleParameter("useLambda")
   force.addGlobalParameter("lambda", lambda_value)
   force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
   force.setCutoffDistance(cutoff) 
   for particle_index in range(nparticles):
      system.addParticle(mass)
      if (particle_index == 0):
         # Add alchemically-modified particle.
         force.addParticle([sigma, epsilon, 1])
      else:
         # Add normal particle.
         force.addParticle([sigma, epsilon, 0])
   system.addForce(force)

   # Add barostat.
   #barostat = MonteCarloBarostat(pressure, temperature)
   #system.addForce(barostat)

   # Store system.
   alchemical_systems.append(system)

# Create random initial positions.
import numpy.random
positions = Quantity(numpy.random.uniform(high=box_edge/angstroms, size=[nparticles,3]), angstrom)

# =============================================================================
# Run simulations at each alchemical lambda value.
# Reduced potentials of sampled configurations are computed for all alchemical states for use with MBAR analysis.
# =============================================================================

u_kln = numpy.zeros([nlambda, nlambda, nprod_iterations], numpy.float64) # u_kln[k,l,n] is reduced potential of snapshot n from simulation k at alchemical state l
for (lambda_index, lambda_value) in enumerate(lambda_values):
   # Create Integrator and Context.
   integrator = LangevinIntegrator(temperature, collision_rate, timestep)
   context = Context(alchemical_systems[lambda_index], integrator)

   # Initiate from last set of positions.
   context.setPositions(positions)

   # Minimize energy from coordinates.
   print("Lambda %d/%d : minimizing..." % (lambda_index+1, nlambda))
   LocalEnergyMinimizer.minimize(context)

   # Equilibrate.
   print("Lambda %d/%d : equilibrating..." % (lambda_index+1, nlambda))
   integrator.step(nequil_steps)

   # Sample.
   position_history = list() # position_history[i] is the set of positions after iteration i
   for iteration in range(nprod_iterations):
      print("Lambda %d/%d : production iteration %d/%d" % (lambda_index+1, nlambda, iteration+1, nprod_iterations))

      # Run dynamics.
      integrator.step(nprod_steps)

      # Get coordinates.
      state = context.getState(getPositions=True)
      positions = state.getPositions(asNumpy=True)

      # Store positions.
      position_history.append(positions)
      
   # Clean up.
   del context, integrator

   # Compute reduced potentials of all snapshots at all alchemical states for MBAR.
   print("Lambda %d/%d : computing energies at all states..." % (lambda_index+1, nlambda))
   beta = 1.0 / (BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA * temperature) # inverse temperature
   for l in range(nlambda):
      # Set up Context just to evaluate energies.
      integrator = VerletIntegrator(timestep)
      context = Context(alchemical_systems[l], integrator)

      # Compute reduced potentials of all snapshots.
      for n in range(nprod_iterations):
         context.setPositions(position_history[n])
         state = context.getState(getEnergy=True)
         u_kln[lambda_index, l, n] = beta * state.getPotentialEnergy()

      # Clean up.
      del context, integrator

# Check to make sure pymbar is installed.
try:
   from timeseries import subsampleCorrelatedData
   from pymbar import MBAR
except:
   raise ImportError("pymbar [https://simtk.org/home/pymbar] must be installed to complete analysis of free energies.")
   
# =============================================================================
# Subsample correlated samples to generate uncorrelated subsample.
# =============================================================================

print("Subsampling data to remove correlation...")
K = nlambda # number of states
N_k = nprod_iterations*numpy.ones([K], numpy.int32) # N_k[k] is the number of uncorrelated samples at state k
u_kln_subsampled = numpy.zeros([K,K,nprod_iterations], numpy.float64) # subsampled data
for k in range(K):
   # Get indices of uncorrelated samples.
   indices = subsampleCorrelatedData(u_kln[k,k,:])
   # Store only uncorrelated data.
   N_k[k] = len(indices)
   for l in range(K):
      u_kln_subsampled[k,l,0:len(indices)] = u_kln[k,l,indices]
print("Number of uncorrelated samples per state:")
print(N_k)

# =============================================================================
# Analyze with MBAR to compute free energy differences and statistical errors.
# =============================================================================

print("Analyzing with MBAR...")
mbar = MBAR(u_kln_subsampled, N_k)
Deltaf_ij, dDeltaf_ij = mbar.getFreeEnergyDifferences()
print("Free energy differences (in kT)")
print(Deltaf_ij)
print("Statistical errors (in kT)")
print(dDeltaf_ij)

# =============================================================================
# Report result.
# =============================================================================

print("Free energy of inserting argon particle: %.3f +- %.3f kT" % (Deltaf_ij[0,K-1], dDeltaf_ij[0,K-1]))

