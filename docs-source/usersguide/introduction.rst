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
http://openmm.org.   Among other things there is a discussion forum,
a wiki, and videos of lectures on using OpenMM.


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

OpenMM software and all related activities, such as this manual, are funded by
the Simbios National Center for Biomedical Computing through the National
Institutes of Health Roadmap for Medical Research, Grant U54 GM072970, and by
National Institutes of Health grant R01-GM062868.

