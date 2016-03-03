## OpenMM: A High Performance Molecular Dynamics Library

Introduction
------------

[OpenMM](http://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. It
provides a combination of extreme flexibility (through custom forces and integrators), openness, and high performance (especially on recent GPUs) that make it truly unique among simulation codes.  This release differs from the 
Main OpenMM release as it contains a softcore VDW for lambda based free energy 
simulations. In addition, we have released an updated version of tinker dynamic_omm.x for
runing simulations using tinker.To install this, move the contents of the online folder http://biomol.bme.utexas.edu/~mh43854/openmm/tinker
 into the source directory of a tinker installation( available at http://dasher.wustl.edu/tinker/).
 Modify the makefile so that all directories point to the appropriate location,
 and then make. documentation specific to this release is available at http://biomol.bme.utexas.edu/wiki/index.php/Research_ex:Openmm


Getting Help
------------

Need Help? Check out the [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).

[C++ API Reference](http://docs.openmm.org/6.3.0/api-c++/namespaceOpenMM.html)

[Python API Reference](http://docs.openmm.org/6.3.0/api-python/annotated.html)

Badges
------
* Travis CI `linux` and `osx` integration tests:
  * GitHub master [![Build Status](https://travis-ci.org/pandegroup/openmm.png?branch=master)](https://travis-ci.org/pandegroup/openmm)
  * `openmm-dev` recipe [![Build Status](https://travis-ci.org/omnia-md/conda-dev-recipes.svg?branch=master)](https://travis-ci.org/omnia-md/conda-dev-recipes)
* Anaconda Cloud `openmm` conda release: [![Binstar `openmm` conda release](https://anaconda.org/omnia/openmm/badges/version.svg)](https://anaconda.org/omnia/openmm)
* Anaconda Cloud `openmm-dev` conda package: [![Binstar `openmm-dev` conda package](https://anaconda.org/omnia/openmm-dev/badges/version.svg)](https://anaconda.org/omnia/openmm-dev)
