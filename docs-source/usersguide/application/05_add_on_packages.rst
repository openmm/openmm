.. default-domain:: py

.. _add-on-packages:

Add-On Packages
###############

In addition to the main OpenMM package, there are several other packages you can
install to add extra features to OpenMM.  In this chapter we will go over some
of the most important ones, describe what they do, and provide references for
where to get more information.

OpenMM-Torch
************

The OpenMM-Torch package provides an interface to the PyTorch machine learning
framework.  It lets you define new types of forces through PyTorch code.  To
install it, open a command-line terminal and type
::

    conda install -c conda-forge openmm-torch

It provides a new class called :class:`TorchForce` that you can add to a System.
The interaction is defined by a PyTorch module that you write using ordinary
Python code, then compile to TorchScript.  The module should have a :meth:`forward`
method that takes the particle positions as input and returns the potential
energy.  Forces are computed through automatic differentiation.

Here is a minimal example of code that applies forces to particles using
PyTorch.  In this case it is a harmonic potential attracting every particle to
the origin.  It could just as easily be any other function of the positions.

.. samepage::
    ::

        import torch
        from openmmtorch import TorchForce

        class ForceModule(torch.nn.Module):
            def forward(self, positions):
                return torch.sum(positions**2)

        module = torch.jit.script(ForceModule())
        force = TorchForce(module)
        system.addForce(force)

    .. caption::

        :autonumber:`Example,OpenMM-Torch example`

OpenMM-Torch has additional features beyond what is shown in this example.  For
example, it supports forces that use periodic boundary conditions, or modules
that directly compute forces rather than relying on automatic differentiation.
For more information, see the OpenMM-Torch_ website.

.. _OpenMM-Torch: https://github.com/openmm/openmm-torch

OpenMM-ML
*********

OpenMM-ML is a higher level utility for modeling molecules with machine learning
potentials.  It provides a selection of standard, pre-trained potential functions
that can be used much like force fields.  To install it, open a command-line
terminal and type
::

    conda install -c conda-forge openmm-ml

It provides a class called :class:`Potential` representing a machine learning
potential function.  To use it, specify the name of the potential function to
use and call :meth:`createSystem` on it to create a System from a Topology.

.. samepage::
    ::

        from openmmml import MLPotential

        potential = MLPotential('ani2x')
        system = potential.createSystem(topology)

    .. caption::

        :autonumber:`Example,OpenMM-ML example`

OpenMM-ML can also create mixed systems that are partly modeled with a machine
learning potential and partly with a conventional force field.  For more
information, see the OpenMM-ML_ website.

.. _OpenMM-ML: https://github.com/openmm/openmm-ml

OpenMMTools
***********

OpenMMTools is a collection of Python utilities providing many useful functions.
Examples include additional integrators, a framework for Markov chain Monte Carlo
simulations, enhanced sampling methods such as replica exchange, and tools for
alchemical simulations to calculate free energies.  To install it, open a
command-line terminal and type
::

    conda install -c conda-forge openmmtools

For more information,  see the OpenMMTools_ website.

.. _OpenMMTools: https://github.com/choderalab/openmmtools

openmmforcefields
*****************

The openmmforcefields package is described in section :numref:`small-molecule-parameters`,
where we show how to use it to parametrize small molecules.  In addition to
its support for generic small molecule force fields like OpenFF and GAFF, it
also provides a larger selection of Amber and CHARMM force fields than what is
bundled with OpenMM.  See the openmmforcefields_ website for more information.

.. _openmmforcefields: http://github.com/openmm/openmmforcefields

OpenMM-PLUMED
*************

PLUMED is a popular tool for computing collective variables of many sorts.  The
OpenMM-PLUMED package allows you to use PLUMED as a force in a simulation.  To
install it, open a command-line terminal and type
::

    conda install -c conda-forge openmm-plumed

PLUMED uses text-based control scripts to define collective variables.  This
package provides a new force class called :class:`PlumedForce` that takes a
PLUMED script and returns the value of the collective variable it defines as its
energy.  For example, this code uses PLUMED to perform metadynamics based on the
distance between particles 0 and 9.


.. samepage::
    ::

        from openmmplumed import PlumedForce

        script = """
        d: DISTANCE ATOMS=1,10
        METAD ARG=d SIGMA=0.2 HEIGHT=0.3 PACE=500"""
        force = PlumedForce(script)
        system.addForce(force)

    .. caption::

        :autonumber:`Example,OpenMM-PLUMED example`

Be aware that PLUMED numbers atoms starting at 1, unlike OpenMM which numbers
them starting at 0.  That is why the script lists the atoms as 1 and 10 rather
than 0 and 9.

For more information, see the websites for the OpenMM-PLUMED_ package and the
PLUMED_ library.

.. _OpenMM-PLUMED: https://github.com/openmm/openmm-plumed
.. _PLUMED: https://www.plumed.org
