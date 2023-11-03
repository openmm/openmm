=============
Extra classes
=============

Tabulated functions
===================

These classes use table of values to define a mathematical function and can be
used by various :ref:`custom forces <custom-forces>`.
The :ref:`OpenMM::TabulatedFunction` class is an abstract class that the other classes
extend.

.. toctree::
    :maxdepth: 2

    generated/TabulatedFunction

    generated/Continuous1DFunction
    generated/Continuous2DFunction
    generated/Continuous3DFunction
    generated/Discrete1DFunction
    generated/Discrete2DFunction
    generated/Discrete3DFunction

Virtual Sites
=============

A virtual site is a particle whose position is computed directly from the
positions of other particles. The :ref:`OpenMM::VirtualSite` class is an abstract
class that the other classes extend.

.. toctree::
    :maxdepth: 2

    generated/VirtualSite

    generated/LocalCoordinatesSite
    generated/OutOfPlaneSite
    generated/ThreeParticleAverageSite
    generated/TwoParticleAverageSite

Serialization
=============

These classes are used to serialize other objects, allowing them to be stored on
disk.

.. toctree::
    :maxdepth: 2

    generated/SerializationNode
    generated/SerializationProxy
    generated/XmlSerializer

Other classes
=============

These classes don't fit neatly into the other categories, but that is not to say
that they aren't important!

.. toctree::
    :maxdepth: 2

    generated/LocalEnergyMinimizer
    generated/MinimizationReporter
    generated/NoseHooverChain
    generated/OpenMMException
    generated/Vec3
