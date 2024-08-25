.. _compiling-openmm-from-source-code:

Compiling OpenMM from Source Code
#################################

This chapter describes the procedure for building and installing OpenMM from
source code.  In most cases, it is best to use the pre-built versions installed
with conda.  Sometimes you might need to build from source, though, such as if
you want to modify OpenMM, or if conda does not provide packages compatible with
your environment.

We first describe how to build on Linux or Mac.  We then describe how to build
on Windows, where the process is slightly different.

.. _compiling-openmm-from-source-linux:

Compiling on Linux and Mac
**************************

Prerequisites
=============

Before building OpenMM from source, you will need certain tools.

C++ compiler
------------

You must have a C++ compiler installed before attempting to build OpenMM from
source.  All recent versions of clang or gcc should work correctly.  On Linux,
you can install the compiler with your system's standard package manager (such
as apt or yum).  On MacOS, you can get a C++ compiler by installing the Xcode
developer tools from the App Store.  Alternatively you can use a package manager
such as Homebrew to install clang or gcc.

Python
------

You will need a 64-bit Python 3.x environment.  We recommend using Miniconda
(https://docs.conda.io/en/latest/miniconda.html), which includes the conda
package manager.

OpenMM Source Code
------------------

You will also need the OpenMM source code.  To download it:

#. Browse to https://github.com/openmm/openmm/releases.
#. Find the latest release and click the link to download the source as either
   a .zip or .tar.gz file.
#. Unpack the file.  Note the location where you unpacked the OpenMM source code.

Alternatively, if you want the most recent development version of the code rather
than the version corresponding to a particular release, you can get it from
https://github.com/openmm/openmm.  Be aware that the development code is constantly
changing, may contain bugs, and should never be used for production work.  If
you want a stable, well tested version of OpenMM, you should download the source
code for the latest release as described above.

CUDA or OpenCL Support
----------------------

If you want to compile OpenMM with support for running on GPUs, you will need
CUDA and/or OpenCL.  MacOS comes with OpenCL built in, so nothing else needs to
be installed.  For Linux, you need an appropriate SDK.

The easiest way is to install the most recent CUDA Toolkit from https://developer.nvidia.com/cuda-downloads.
It includes the headers and libraries needed to compile both CUDA and OpenCL
applications.  In addition, it has runtime libraries that are needed for running
CUDA applications.  The runtime components for OpenCL applications are included
with the GPU drivers from NVIDIA, AMD, and Intel, so make sure you have an
up-to-date driver.

Other Required Software
-----------------------

Several other tools are required to build OpenMM.  The easiest way to install
them is with conda.  The following command will install everything needed to
build OpenMM.
::

    conda install -c conda-forge cmake make cython swig doxygen numpy

Step 1: Configure with CMake
============================

First, create a directory in which to build OpenMM.  Starting from the root
level of the OpenMM source tree (the directory containing the top CMakeLists.txt
file), execute the following commands:

.. code-block:: none

    mkdir build
    cd build
    ccmake ..

That is not a typo.  :code:`ccmake` has two c’s.  CCMake is the visual CMake
configuration tool.  Press :code:`c` within the CCMake interface to load the
OpenMM build scripts and begin configuring CMake.

There are several variables that can be adjusted in the CMake interface:

* Set the variable CMAKE_INSTALL_PREFIX to the location where you want to
  install OpenMM. If you are using conda environments this variable should point to
  the full path of the root directory of your environment.
* Set the variable PYTHON_EXECUTABLE to the Python interpreter you plan to use
  OpenMM with.  Usually this will be detected automatically.
* There are lots of options starting with OPENMM_BUILD that control
  whether to build particular features of OpenMM, such as plugins, API wrappers,
  and documentation.
* Usually the OpenCL library and headers will be detected automatically.  If for
  any reason CMake is unable to find them, set OPENCL_INCLUDE_DIR to point to
  the directory containing the headers (usually /usr/local/cuda/include on Linux)
  and OPENCL_LIBRARY to point to the library (usually /usr/local/cuda/lib64/libOpenCL.so
  on Linux).

Configure (press “c”) again.  Adjust any variables that cause an error.

Continue to configure (press “c”) until no starred CMake variables are
displayed, then press “g” to generate the makefiles for building the project.

Step 2: Build
=============

Build OpenMM with the command::

    make

.. _test-your-build:

Step 3: Test your build
=======================

This step is optional but recommended. Tests can take up to several minutes depending on your
hardware configuration.

It is recommended that you make sure your local build of OpenMM works before trying
to install.

After OpenMM has been built, you should run the unit tests to make sure it
works.  Enter the command::

    make test

You should see a series of test results like this:

.. code-block:: none

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

:code:`Passed` is good.  :code:`FAILED` is bad.  If any tests fail, you
can run them individually to get more detailed error information.  For example,
::

    ./TestReferenceLangevinIntegrator

Note that some tests are stochastic, and therefore are expected to fail a small
fraction of the time.  These tests will say so in the error message:

.. code-block:: none

    exception: Assertion failure at TestReferenceLangevinIntegrator.cpp:129.  Expected 9.97741,
        found 10.7884 (This test is stochastic and may occasionally fail)

Step 3: Install
===============
Install your local build of OpenMM using the following command::

    make install

If you are installing to a system directory, such as /usr/local/openmm/, you will
need admin capabilities to install, in this case use::

    sudo make install

Step 3: Install the Python API
==============================

Build and install the Python API with the command::

    make PythonInstall

If you are installing into the system Python, such as /usr/bin/python, you will
need to type::

    sudo make PythonInstall

You can test the Python API installation using::

    python -m openmm.testInstallation

Congratulations! You have successfully built and installed OpenMM from source.

Compiling on Windows
********************

Prerequisites
=============

Before building OpenMM from source, you will need certain tools.

C++ compiler
------------

On Windows systems, use the C++ compiler in Visual Studio 2017 or later.  You
can download a free version of Visual Studio from https://visualstudio.microsoft.com.

Python
------

You will need a 64-bit Python 3.x environment.  We recommend using Miniconda
(https://docs.conda.io/en/latest/miniconda.html), which includes the conda
package manager.

CMake
-----

CMake is the build system used for OpenMM.  You must install CMake version 3.17
or higher before attempting to build OpenMM from source.  You can get CMake from
http://www.cmake.org/.

OpenMM Source Code
------------------

You will also need the OpenMM source code.  To download it:

#. Browse to https://github.com/openmm/openmm/releases.
#. Find the latest release and click the link to download the source as either
   a .zip or .tar.gz file.
#. Unpack the file.  Note the location where you unpacked the OpenMM source code.

Alternatively, if you want the most recent development version of the code rather
than the version corresponding to a particular release, you can get it from
https://github.com/openmm/openmm.  Be aware that the development code is constantly
changing, may contain bugs, and should never be used for production work.  If
you want a stable, well tested version of OpenMM, you should download the source
code for the latest release as described above.

CUDA or OpenCL Support
----------------------

If you want to compile OpenMM with support for running on GPUs, you will need
CUDA and/or OpenCL.  Install the most recent CUDA Toolkit from https://developer.nvidia.com/cuda-downloads.
It includes the headers and libraries needed to compile both CUDA and OpenCL
applications.  In addition, it has runtime libraries that are needed for running
CUDA applications.  The runtime components for OpenCL applications are included
with the GPU drivers from NVIDIA, AMD, and Intel, so make sure you have an
up-to-date driver.

Other Required Software
-----------------------

Several other tools are required to build OpenMM.  The easiest way to install
them is with conda.  From the Windows Start menu, select "Anaconda Prompt (Miniconda3)".
It will open a command window that is preconfigured for conda.  Enter the
following command to install everything needed to build OpenMM.
::

    conda install -c conda-forge cython swig doxygen numpy

Step 1: Configure with CMake
============================

First, create a directory in which to build OpenMM.  In the "Anaconda Prompt"
window opened above, cd to the root level of the OpenMM source tree (the
directory containing the top CMakeLists.txt file).  Execute the following commands:

.. code-block:: none

    mkdir build
    cd build
    "C:\Program Files\CMake\bin\cmake-gui.exe"

This will launch the CMake GUI configuration tool.  It is critical that you
launch it from the "Anaconda Prompt" window as shown above.  Do *not* launch
it from the Start menu.  If you do, it will not be able to find the tools you
installed with conda.

#. In the box labeled "Where is the source code" browse to the OpenMM source directory
   (containing the top CMakeLists.txt file).
#. In the box labeled "Where to build the binaries" browse to the build
   directory you just created.
#. Click the "Configure" button at the bottom of the CMake window.
#. Select "Visual Studio 16 2019" from the  list of Generators (or whichever
   version you have installed) and click "Finish".

There are several variables that can be adjusted in the CMake interface:

* Set the variable CMAKE_INSTALL_PREFIX to the location where you want to
  install OpenMM.
* Set the variable PYTHON_EXECUTABLE to the Python interpreter you plan to use
  OpenMM with.  Usually this will be detected automatically.
* There are lots of options starting with OPENMM_BUILD that control
  whether to build particular features of OpenMM, such as plugins, API wrappers,
  and documentation.
* Usually the OpenCL library and headers will be detected automatically.  If for
  any reason CMake is unable to find them, set OPENCL_INCLUDE_DIR to point to
  the directory containing the headers (usually
  "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/include", except
  with the correct version number for the toolkit you installed) and
  OPENCL_LIBRARY to point to the library (usually "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/lib/x64/OpenCL.lib").

Press "Configure" again.  Adjust any variables that cause an error.

Continue to press "Configure" until no red CMake variables are displayed, then
press "Generate" to create the Visual Studio project files for building OpenMM.

Step 2: Build and Install
=========================

#. Open the file :file:`OpenMM.sln` in your build directory in Visual Studio.
   Note that this file will appear as just :file:`OpenMM` if you have configured
   Explorer to hide file name extensions.
#. Set the configuration type to "Release" (not "Debug") in the toolbar.
#. From the Build menu, select "Build Solution".  This takes some time.
#. In the Solution Explorer, right-click on "INSTALL" and select "Build".

Step 3: Install the Python API
==============================

In the Solution Explorer, right-click on "PythonInstall" and select "Build".

Step 4: Test your build
=======================

After OpenMM has been built, you should run the unit tests to make sure it
works.  In the Solution Explorer, right-click on "RUN_TESTS" and select "Build".
You should see a series of test results like this:

.. code-block:: none

            Start   1: TestReferenceAndersenThermostat
      1/317 Test   #1: TestReferenceAndersenThermostat .............. Passed  0.26 sec
            Start   2: TestReferenceBrownianIntegrator
      2/317 Test   #2: TestReferenceBrownianIntegrator .............. Passed  0.13 sec
            Start   3: TestReferenceCheckpoints
      3/317 Test   #3: TestReferenceCheckpoints ..................... Passed  0.02 sec
      ... <many other tests> ...

:code:`Passed` is good.  :code:`FAILED` is bad.  If any tests fail, you
can run them individually to get more detailed error information.  Right-click
on a test in the Solution Explorer and select "Debug > Start New Instance".

Note that some tests are stochastic, and therefore are expected to fail a small
fraction of the time.  These tests will say so in the error message:

.. code-block:: none

    exception: Assertion failure at TestReferenceLangevinIntegrator.cpp:129.  Expected 9.97741,
        found 10.7884 (This test is stochastic and may occasionally fail)

Congratulations! You have successfully built and installed OpenMM from source.


Building the Documentation (Optional)
*************************************

The documentation that you're currently reading, as well as the Developer Guide and API
documentation, can be built through CMake.  To do that, you need to install a few
additional tools.  The easiest way is to use :code:`conda` to install them into
your Python environment.  The following command installs everything needed to
build documentation in HTML format.
::

    conda install -c conda-forge sphinx sphinxcontrib-bibtex breathe jinja2

To build documentation in PDF format, you also must have a functional LaTeX
installation.  It can be obtained from https://www.latex-project.org/get.

If you want to build documentation, make sure that OPENMM_GENERATE_API_DOCS is
set to ON when configuring the build in CMake.

To build the documentation, use the following build targets.

* :code:`sphinxhtml`: Build the User Guide and Developer Guide in HTML format.

* :code:`sphinxpdf`: Build the User Guide and Developer Guide in PDF format.

* :code:`C++ApiDocs`: Build the C++ API documentation.

* :code:`PythonApiDocs`: Build the Python API documentation.  This target
  requires that you have already built the :code:`install` target, such as with
  :code:`make install`.

On Linux or Mac, build a target using the :code:`make` command.  For example,
::

    make sphinxhtml

On Windows, right-click on the target in the Solution Explorer and select "Build".

After building the documentation, build the :code:`install` target to install
the documentation into the installation directory (the one you specified with
CMAKE_INSTALL_PREFIX).

Using local build of OpenMM alongside conda tools that depend on it
*******************************************************************

A common case is to have a local build of OpenMM in the same environment as other tools
that depend on it. This can be achieved by forcing a remove of OpenMM when you install
your tools using conda.

We will use :code:`openmmtools` as an example here, but it can be replaced with any
other software package that requires OpenMM.

Step 1: Install your tools as usual
===================================

Install your tools using conda as you commonly do, for example using::

    conda install -c conda-forge  openmmtools

This will pull the conda-forge package of :code:`openmm` which we don't want since we want
to use our local build.

Step 2: Remove conda openmm package
===================================

To remove the openmm package that was installed in the previous step, we can use::

    conda remove --force openmm

This will remove the :code:`openmm` package without changing or removing dependencies.

Step 3: Install local build of openmm
=====================================

Now we just install our local build of :code:`openmm` as instructed in
:ref:`compiling-openmm-from-source-linux`
