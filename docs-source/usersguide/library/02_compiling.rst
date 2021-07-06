.. _compiling-openmm-from-source-code:

Compiling OpenMM from Source Code
#################################

This chapter describes the procedure for building and installing OpenMM
libraries from source code.  It is recommended that you use binary OpenMM
libraries, if possible.  If there are not suitable binary libraries for your
system, consider building OpenMM from source code by following these
instructions.

Prerequisites
*************

Before building OpenMM from source, you will need the following:

* A C++ compiler

* CMake

* OpenMM source code


See the sections below for specific instructions for the different platforms.

Get a C++ compiler
==================

You must have a C++ compiler installed before attempting to build OpenMM from
source.

Mac and Linux: clang or gcc
---------------------------

Use clang or gcc on Mac/Linux.  OpenMM should compile correctly with all recent
versions of these compilers.  We recommend clang since it produces faster code,
especially when using the CPU platform.  If you do not already have a compiler
installed, you will need to download and install it.  On Mac OS X, this means
downloading the Xcode Tools from the App Store.

Windows: Visual Studio
----------------------

On Windows systems, use the C++ compiler in Visual Studio 2015 or later.  You
can download a free version of the Visual Studio C++ build tools from
https://visualstudio.microsoft.com.  If you plan to use OpenMM from
Python, it is critical that both OpenMM and Python be compiled with the same
version of Visual Studio.

Install CMake
=============

CMake is the build system used for OpenMM.  You must install CMake version 3.1
or higher before attempting to build OpenMM from source.  You can get CMake from
http://www.cmake.org/.  If you choose to build CMake from source on Linux, make
sure you have the curses library installed beforehand, so that you will be able
to build the CCMake visual CMake tool.

Get the OpenMM source code
==========================

You will also need the OpenMM source code before building OpenMM from source.
To download and unpack OpenMM source code:

#. Browse to https://simtk.org/home/openmm.
#. Click the "Downloads" link in the navigation bar on the left side.
#. Download OpenMM<Version>-Source.zip, choosing the latest version.
#. Unpack the zip file.  Note the location where you unpacked the OpenMM source
   code.

Alternatively, if you want the most recent development version of the code rather
than the version corresponding to a particular release, you can get it from
https://github.com/pandegroup/openmm.  Be aware that the development code is constantly
changing, may contain bugs, and should never be used for production work.  If
you want a stable, well tested version of OpenMM, you should download the source
code for the latest release as described above.

Other Required Software
=======================

There are several other pieces of software you must install to compile certain
parts of OpenMM.  Which of these you need depends on the options you select in
CMake.

* For compiling the CUDA Platform, you need:

   * CUDA (See Chapter :numref:`installing-openmm` for installation instructions.)

* For compiling the OpenCL Platform, you need:

   * OpenCL (See Chapter :numref:`installing-openmm` for installation instructions.)

* For compiling C and Fortran API wrappers, you need:

   * Python 2.7 or later (http://www.python.org)
   * Doxygen (http://www.doxygen.org)
   * A Fortran compiler

* For compiling the Python API wrappers, you need:

   * Python 2.7 or later (http://www.python.org)
   * SWIG (http://www.swig.org)
   * Doxygen (http://www.doxygen.org)
   * Cython (https://cython.org)

* For compiling the CPU platform, you need:

   * FFTW, single precision multithreaded version (http://www.fftw.org)

* To generate API documentation, you need:

   * Doxygen (http://www.doxygen.org)



Step 1: Configure with CMake
****************************


Build and source directories
============================

First, create a directory in which to build OpenMM.  A good name for this
directory is build_openmm.  We will refer to this as the “build_openmm
directory” in the instructions below.  This directory will contain the temporary
files used by the OpenMM CMake build system.  Do not create this build directory
within the OpenMM source code directory.  This is what is called an “out of
source” build, because the build files will not be mixed with the source files.

Also note the location of the OpenMM source directory (i.e., where you unpacked
the source code zip file).  It should contain a file called CMakeLists.txt.
This directory is what we will call the “OpenMM source directory” in the
following instructions.

Starting CMake
==============

Configuration is the first step of the CMake build process.  In the
configuration step, the values of important build variables will be established.

Mac and Linux
-------------

On Mac and Linux machines, type the following two lines:
::

    cd build_openmm
    ccmake -i <path to OpenMM src directory>

That is not a typo.  :code:`ccmake` has two c’s.  CCMake is the visual CMake
configuration tool.         Press “\ :code:`c`\ ” within the CCMake interface to
configure CMake.  Follow the instructions in the “All Platforms” section below.

Windows
-------

On Windows, perform the following steps:

#. Click Start->All Programs->CMake 3.1->CMake
#. In the box labeled "Where is the source code:" browse to OpenMM src directory
   (containing top CMakeLists.txt)
#. In the box labeled "Where to build the binaries" browse to your build_openmm
   directory.
#. Click the "Configure" button at the bottom of the CMake screen.
#. Select "Visual Studio 10 2010" from the  list of Generators (or whichever
   version you have installed)
#. Follow the instructions in the “All Platforms” section below.


All platforms
-------------

There are several variables that can be adjusted in the CMake interface:

* If you intend to use CUDA (NVIDIA) or OpenCL acceleration, set the variable
  OPENMM_BUILD_CUDA_LIB or OPENMM_BUILD_OPENCL_LIB, respectively, to ON.  Before
  doing so, be certain that you have installed and tested the drivers for the
  platform you have selected (see Chapter :numref:`installing-openmm` for information on
  installing GPU software).
* There are lots of other options starting with OPENMM_BUILD that control
  whether to build particular features of OpenMM, such as plugins, API wrappers,
  and documentation.
* Set the variable CMAKE_INSTALL_PREFIX to the location where you want to
  install OpenMM.
* Set the variable PYTHON_EXECUTABLE to the Python interpreter you plan to use
  OpenMM with.


Configure (press “c”) again.  Adjust any variables that cause an
error.

Continue to configure (press “c”) until no starred/red CMake
variables are displayed.  Congratulations, you have completed the configuration
step.

Step 2: Generate Build Files with CMake
***************************************

Once the configuration is done, the next step is generation.  The generate
“g” or “OK” or “Generate” option will not be available until
configuration has completely converged.

Windows
=======

* Press the "OK" or “Generate” button to generate Visual Studio project files.

* If CMake does not exit automatically, press the close button in the upper-
  right corner of the CMake title bar to exit.


Mac and Linux
=============

* Press “g” to generate the Makefile.
* If CMake does not exit automatically, press “q” to exit.


That’s it!  Generation is the easy part.  Now it’s time to build.

Step 3: Build OpenMM
********************


Windows
=======

#. Open the file OpenMM.sln in your openmm_build directory in Visual Studio.
#. Set the configuration type to "Release" (not "Debug") in the toolbar.
#. From the Build menu, click Build->Build Solution
#. The OpenMM libraries and test programs will be created.  This takes some
   time.
#. The test program TestCudaRandom might not build on Windows.  This is OK.


Mac and Linux
=============

* Type :code:`make` in the openmm_build directory.

The OpenMM libraries and test programs will be created.  This takes some time.


Step 4: Install OpenMM
**********************


Windows
=======

In the Solution Explorer Panel, far-click/right-click INSTALL->build.

Mac and Linux
=============

Type:
::

    make install

If you are installing to a system area, such as /usr/local/openmm/, you will
need to type:
::

    sudo make install

Step 5: Install the Python API
******************************


Windows
=======

In the Solution Explorer Panel, right-click PythonInstall->build.

Mac and Linux
=============

Type:
::

    make PythonInstall

If you are installing into the system Python, such as /usr/bin/python, you will
need to type:
::

    sudo make PythonInstall

.. _test-your-build:

Step 6: Test your build
***********************

After OpenMM has been built, you should run the unit tests to make sure it
works.

Windows
=======

In Visual Studio, far-click/right-click RUN_TESTS in the Solution Explorer
Panel.  Select RUN_TESTS->build to begin testing.  Ignore any failures for
TestCudaRandom.

Mac and Linux
=============

Type:
::

    make test

You should see a series of test results like this:
::

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

:code:`Passed` is good.  :code:`FAILED` is bad.  If any tests fail, you
can run them individually to get more detailed error information.  Note that
some tests are stochastic, and therefore are expected to fail a small fraction
of the time.  These tests will say so in the error message:
::

    ./TestReferenceLangevinIntegrator

    exception: Assertion failure at TestReferenceLangevinIntegrator.cpp:129.  Expected 9.97741,
        found 10.7884 (This test is stochastic and may occasionally fail)

Congratulations! You successfully have built and installed OpenMM from source.

Building the Documentation (Optional)
*************************************

The documentation that you're currently reading, as well as the developer guide and API
documentation can be built through CMake by setting the OpenMM option :code:`OPENMM_GENERATE_API_DOCS=ON`.

User Guide and Developer Guide
==============================

Generating the user guide and developer guide requires the following dependencies

* Sphinx (http://sphinx-doc.org/)

* sphinxcontrib-bibtex (https://pypi.python.org/pypi/sphinxcontrib-bibtex)

These dependencies may not be available in your system package manager, but should
be installable through Python's ``pip`` package manager. ::

   pip install sphinx sphinxcontrib-bibtex

The developer and user guides can be built either as HTML or a PDFs. Building the
PDF version will also require a functional LaTeX installation.

To build the HTML version of the documentation, type: ::

  make sphinxhtml

To build the PDF version of the documentation, type: ::

  make sphinxpdf


Python and C++ API Documentation
================================

The following dependencies are required to build the Python and C++ API documentation.

* Sphinx (http://sphinx-doc.org/)

* sphinxcontrib-lunrsearch (https://pypi.python.org/pypi/sphinxcontrib-lunrsearch)

* sphinxcontrib-autodoc_doxygen (https://pypi.python.org/pypi/sphinxcontrib-autodoc_doxygen)


These dependencies may not be available in your system package manager, but should
be installable through Python's ``pip`` package manager. ::

   pip install sphinx sphinxcontrib-lunrsearch sphinxcontrib-autodoc_doxygen

To build the C++ API documentation, type: ::

  make C++ApiDocs

To build the Python API documentation, type: ::

  make PythonApiDocs

