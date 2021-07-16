=============
Extra classes
=============

Tabulated functions
===================

These classes use table of values to define a mathematical function and can be
used by various `custom forces <my-reference-label>`_.
The :doc:`generated/OpenMM::TabulatedFunction` class is an abstract class that the other classes
extend.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::TabulatedFunction

    generated/OpenMM::Continuous1DFunction
    generated/OpenMM::Continuous2DFunction
    generated/OpenMM::Continuous3DFunction
    generated/OpenMM::Discrete1DFunction
    generated/OpenMM::Discrete2DFunction
    generated/OpenMM::Discrete3DFunction

Virtual Sites
=============

A virtual site is a particle whose position is computed directly from the
positions of other particles. The :doc:`generated/OpenMM::VirtualSite` class is an abstract
class that the other classes extend.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::VirtualSite

    generated/OpenMM::LocalCoordinatesSite
    generated/OpenMM::OutOfPlaneSite
    generated/OpenMM::ThreeParticleAverageSite
    generated/OpenMM::TwoParticleAverageSite

Serialization
=============

These classes are used to serialize other objects, allowing them to be stored on
disk.

.. toctree::
    :maxdepth: 2

    generated/OpenMM::SerializationNode
    generated/OpenMM::SerializationProxy
    generated/OpenMM::XmlSerializer

Other classes
=============

These classes don't fit neatly into the other categories, but that is not to say
that they aren't important!

.. toctree::
    :maxdepth: 2

    generated/OpenMM::LocalEnergyMinimizer
    generated/OpenMM::NoseHooverChain
    generated/OpenMM::OpenMMException
    generated/OpenMM::Vec3
