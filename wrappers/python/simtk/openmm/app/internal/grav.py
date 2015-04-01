"""
Module for adding gravitational forces to a created system
"""
from simtk import openmm as mm
from simtk import unit as u

import warnings

def addGravity(system, elevation=None):
    """ Adds a new force to model gravitational interaction between particles in
    a system as well as, optionally, the gravitational effect of earth on each
    particle

    Parameters
    ----------
    system : openmm.System
        The system to add a gravitational force to
    elevation : length Quantity, optional
        The elevation above Earth's surface of the *bottom* of the container for
        this system that this simulation is to take place. If None, Earth's
        gravity is ignored (this is NOT recommended for optimal accuracy). We
        assume the Z-dimension is up
    """
    # Define the gravitational constant. Taken not from unreliable sources like
    # the CRC or NIST, but rather from the most reliable of sources: Wikipedia
    G = 6.673848e-11 * u.newtons * u.meters**2 / u.kilograms**2
    # Define the energy function (strictly attractive) of gravity, where each
    # particle is defined only by their mass
    intergrav = mm.CustomNonbondedForce('-G*mass1*mass2/r; G=%g;' %
                        G.value_in_unit_system(u.md_unit_system))
    intergrav.addPerParticleParameter('mass')
    # Add all of the particles' masses to the gravitational force, then add it
    # to the system.
    for i in range(system.getNumParticles()):
        mass = system.getParticleMass()
        intergrav.addParticle((mass,))
    system.addForce(intergrav)
    # Add the force due to gravity of the earth, if requested. It is
    # approximately constant over a distance of nanometers compared to earth's
    # radius.
    if elevation is None:
        warnings.warn('You are ignoring the force of gravity due to the earth. '
                      'The force field will not be as accurate as it could be')
        return
    if u.is_quantity(elevation):
        elevation = elevation.value_in_unit(u.meters)
    # Don't support simulations below sea level yet
    if elevation <= 0:
        raise ValueError('Simulations below sea level are not yet supported')
    elevation *= u.meters
    # Calculate the acceleration; assume Earth's radius is 6.37101e6 meters,
    # mass is 5.9736e24 kg, and g is given above
    g = G * 5.9736e24*u.kilograms / (6.37101e6*u.meters + elevation)**2
    # Now add the external force
    earthgrav = mm.CustomExternalForce('mass*g*z; g=%g' %
                                       g.value_in_unit_system(u.md_unit_system))
    earthgrav.addPerParticleParameter('mass')
    for i in range(system.getNumParticles()):
        mass = system.getParticleMass()
        earthgrav.addParticle((mass,))
    system.addForce(earthgrav)
