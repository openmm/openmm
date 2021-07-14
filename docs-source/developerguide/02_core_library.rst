.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-core-library:

The Core Library
################

OpenMM is based on a layered architecture, as shown in the following diagram:

.. figure:: ../images/ArchitectureLayers.jpg
   :align: center
   :width: 100%

   :autonumber:`Figure,Architecture Layers`\ : OpenMM architecture

The public API layer consists of the classes you access when using OpenMM in an
application: System; Force and its subclasses; Integrator and its subclasses;
and Context.  These classes define a public interface but do no computation.

The next layer down consists of “implementation” classes that mirror the public
API classes: ContextImpl, ForceImpl, and a subclass of ForceImpl for each
subclass of Force (HarmonicBondForceImpl, NonbondedForceImpl, etc.).  These
objects are created automatically when you create a Context.  They store
information related to a particular simulation, and define methods for
performing calculations.

Note that, whereas a Force is logically “part of” a System, a ForceImpl is
logically “part of” a Context.  (See :autonumref:`Figure,API Relationships`\ .)  If you create many Contexts
for simulating the same System, there is still only one System and only one copy
of each Force in it.  But there will be separate ForceImpls for each Context,
and those ForceImpls store information related to their particular Contexts.


.. figure:: ../images/SystemContextRelationships.jpg
   :align: center

   :autonumber:`Figure,API Relationships`\ : Relationships between public API and implementation layer objects

Also note that there is no “IntegratorImpl” class, because it is not needed.
Integrator is already specific to one Context.  Many Contexts can all simulate
the same System, but each of them must have its own Integrator, so information
specific to one simulation can be stored directly in the Integrator.

The next layer down is the OpenMM Low Level API (OLLA).  The important classes
in this layer are: Platform; Kernel; KernelImpl and its subclasses; and
KernelFactory.  A Kernel is just a reference counted pointer to a KernelImpl;
the real work is done by KernelImpl objects (or more precisely, by instances of
its subclasses).  A KernelFactory creates KernelImpl objects, and a Platform
ties together a set of KernelFactories, as well as defining information that
applies generally to performing computations with that Platform.

All of these classes (except Kernel) are abstract.  A particular Platform
provides concrete subclasses of all of them.  For example, the reference
platform defines a Platform subclass called ReferencePlatform, a KernelFactory
subclass called ReferenceKernelFactory, and a concrete subclass of each abstract
KernelImpl type: ReferenceCalcNonbondedForceKernel extends
CalcNonbondedForceKernel (which in turn extends KernelImpl),
ReferenceIntegrateVerletStepKernel extends IntegrateVerletStepKernel, and so on.

We can understand this better by walking through the entire sequence of events
that takes place when you create a Context.  As an example, suppose you create a
System; add a NonbondedForce to it; create a VerletIntegrator; and then create a
Context for them using the reference Platform.  Here is what happens.

#. The Context constructor creates a ContextImpl.
#. The ContextImpl calls :code:`createImpl()` on each Force in the System,
   which creates an instance of the appropriate ForceImpl subclass.
#. The ContextImpl calls :code:`contextCreated()` on the Platform(), which
   in turn calls :code:`setPlatformData()` on the ContextImpl.  This allows
   Platform-specific information to be stored in a ContextImpl.  Every Platform has
   its own mechanism for storing particle masses, constraint definitions, particle
   positions, and so on.  ContextImpl therefore allows the Platform to create an
   arbitrary block of data and store it where it can be accessed by that Platform’s
   kernels.
#. The ContextImpl  calls :code:`createKernel()` on the Platform several
   times to get instances of various kernels that it needs:
   CalcKineticEnergyKernel, ApplyConstraintsKernel, etc.

   #. For each kernel, the Platform looks up which KernelFactory has been
      registered for that particular kernel.  In this case, it will be a
      ReferenceKernelFactory.
   #. It calls :code:`createKernelImpl()` on the KernelFactory, which
      creates and returns an instance of an appropriate KernelImpl subclass:
      ReferenceCalcKineticEnergyKernel, ReferenceApplyConstraintsKernel, etc.

#. The ContextImpl loops over all of its ForceImpls and calls
   :code:`initialize()` on each one.

   #. Each ForceImpl asks the Platform to create whatever kernels it needs.  In
      this example, NonbondedForceImpl will request a CalcNonbondedForceKernel, and
      get back a ReferenceCalcNonbondedForceKernel.

#. The ContextImpl calls :code:`initialize()` on the Integrator which, like
   the other objects, requests kernels from the Platform.  In this example,
   VerletIntegrator requests an IntegrateVerletStepKernel and gets back a
   ReferenceIntegrateVerletStepKernel.


At this point, the Context is fully initialized and ready for doing computation.
Reference implementations of various KernelImpls have been created, but they are
always referenced through abstract superclasses.  Similarly, data structures
specific to the reference Platform have been created and stored in the
ContextImpl, but the format and content of these structures is opaque to the
ContextImpl.  Whenever it needs to access them (for example, to get or set
particle positions), it does so through a kernel (UpdateStateDataKernel in this
case).

Now suppose that you call :code:`step()` on the VerletIntegrator.  Here is
what happens to execute each time step.

#. The VerletIntegrator calls :code:`updateContextState()` on the
   ContextImpl.  This gives each Force an opportunity to modify the state of the
   Context at the start of each time step.

   #. The ContextImpl loops over its ForceImpls and calls
      :code:`updateContextState()` on each one.  In this case, our only ForceImpl is
      a NonbondedForceImpl, which returns without doing anything.  On the other hand,
      if we had an AndersenThermostat in our System, its ForceImpl would invoke a
      kernel to modify particle velocities.

#. The VerletIntegrator calls :code:`calcForcesAndEnergy()` on the
   ContextImpl to request that the forces be computed.

   #. The ContextImpl calls :code:`beginComputation()` on its
      CalcForcesAndEnergyKernel.  This initializes all the forces to zero and does any
      other initialization the Platform requires before forces can be computed.  For
      example, some Platforms construct their nonbonded neighbor lists at this point.
   #. The ContextImpl loops over its ForceImpls and calls
      :code:`calcForcesAndEnergy()` on each one.  In this case, we have a
      NonbondedForceImpl which invokes its CalcNonbondedForceKernel to compute forces.
   #. Finally, the ContextImpl calls :code:`finishComputation()` on its
      CalcForcesAndEnergyKernel.  This does any additional work needed to determine
      the final forces, such as summing the values from intermediate buffers.

#. Finally, the VerletIntegrator invokes its IntegrateVerletStepKernel.  This
   takes the forces, positions, and velocities that are stored in a Platform-
   specific format in the ContextImpl, uses them to compute new positions and
   velocities, and stores them in the ContextImpl.

