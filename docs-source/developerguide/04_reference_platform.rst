.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++


.. _the-reference-platform:

The Reference Platform
######################

The reference Platform is written with simplicity and clarity in mind, not
performance.  (It is still not always as simple or clear as one might hope, but
that is the goal.)  When implementing a new feature, it is recommended to create
the reference implementation first, then use that as a model for the versions in
other Platforms.

When using the reference Platform, the “platform-specific data” stored in
ContextImpl is of type ReferencePlatform::PlatformData, which is declared in
ReferencePlatform.h.  It has fields for storing positions, velocities, box
vectors, and other types of data.

The PlatformData’s vector of forces contains one element for each particle.  At
the start of each force evaluation, all elements of it are set to zero.  Each
Force adds its own contributions to the vector, so that at the end, it contains
the total force acting on each particle.

There are a few additional classes that contain useful static methods.
SimTKOpenMMUtilities has various utility functions, of which the most important
is a random number generator.  ReferenceForce provides methods for calculating
the displacement between two positions, optionally taking periodic boundary
conditions into account.

