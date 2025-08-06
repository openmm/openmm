Introduction
############

OpenMM consists of two parts:

#. A set of libraries that lets programmers easily add molecular simulation
   features to their programs
#. An “application layer” that exposes those features to end users who just want
   to run simulations


This guide is divided into three sections:

* :ref:`Part I <the-openmm-application-layer>`
  describes the application layer.  It is relevant to all users, but especially relevant to people
  who want to use OpenMM as a stand-alone application for running simulations.
* :ref:`Part II <the-openmm-library>`
  describes how to use the OpenMM libraries within your own applications.  It is primarily
  relevant to programmers who want to write simulation applications.
* :ref:`Part III <the-theory-behind-openmm>`
  describes the mathematical theory behind the features found in OpenMM.  It is relevant to all users.


Online Resources
****************

You can find more documentation and other material at our website
https://openmm.org/.  In particular, you may want to consult the
`Python API documentation <https://docs.openmm.org/latest/api-python/>`__ or the
`C++ API documentation <https://docs.openmm.org/latest/api-c++/>`__, read the
`Frequently Asked Questions (FAQ) <https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions>`__,
or browse the `discussion forum on GitHub <https://github.com/openmm/openmm/discussions>`__.


Referencing OpenMM
******************

Any work that uses OpenMM should cite the following publication:

P. Eastman, R. Galvelis, R. P. Peláez, C. R. A. Abreu, S. E. Farr, E. Gallicchio,
A. Gorenko, M. M. Henry, F. Hu, J. Huang, A. Krämer, J. Michel, J. A. Mitchell,
V. S. Pande, J. PGLM Rodrigues, J. Rodriguez-Guerra, A. C. Simmonett, S. Singh,
J. Swails, P. Turner, Y. Wang, I. Zhang, J. D. Chodera, G. De Fabritiis, and
T. E. Markland. "OpenMM 8: Molecular Dynamics Simulation with Machine Learning
Potentials." J. Phys. Chem. B 128(1), pp. 109-116 (2023).

We depend on academic research grants to fund the OpenMM development efforts;
citations of our publication will help demonstrate the value of OpenMM.


Acknowledgments
***************

OpenMM research and development activities are supported by various individuals
and funded from a number of sources: for up-to-date information, consult
https://openmm.org/development.
