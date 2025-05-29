.. default-domain:: py

.. _the-openmm-application-layer-introduction:

Getting Started
###############

Introduction
************

The first thing to understand about the OpenMM “application layer” is that it is
not exactly an application in the traditional sense: there is no program called
“OpenMM” that you run.  Rather, it is a collection of libraries written in the
Python programming language.  Those libraries can easily be chained together to
create Python programs that run simulations.  But don’t worry!  You don’t need
to know anything about Python programming (or programming at all) to use it.
Nearly all molecular simulation applications ask you to write some sort of
“script” that specifies the details of the simulation to run.  With OpenMM, that
script happens to be written in Python.  But it is no harder to write than those
for most other applications, and this guide will teach you everything you need
to know.  There is even a graphical interface that can write the script for you
based on a simple set of options (see Section :numref:`the-script-builder-application`),
so you never need to type a single line of code!

On the other hand, if you don’t mind doing a little programming, this approach
gives you enormous power and flexibility.  Your script has complete access to
the entire OpenMM application programming interface (API), as well as the full
power of the Python language and libraries.  You have complete control over
every detail of the simulation, from defining the molecular system to analyzing
the results.


.. _installing-openmm:

Installing OpenMM
*****************

OpenMM is installed using either the pip or conda package manager.  Pip is
normally included with all Python distributions.  The easiest way to get conda
is to use Miniconda (available from https://docs.conda.io/en/latest/miniconda.html),
a Python distribution with conda preinstalled.

(Another option is to compile OpenMM from source.  This provides more flexibility,
but it is much more work, and there is rarely a need for anyone but advanced users
to compile from source.  Detailed instruction are in Chapter :numref:`compiling-openmm-from-source-code`.)

\1. Begin by installing the most recent 64 bit, Python 3.x version of Miniconda.

\2. (Optional) If you want to run OpenMM on a GPU, make sure you have installed
modern drivers from your vendor.

  * If you have an NVIDIA GPU, download the latest drivers from
    https://www.nvidia.com/Download/index.aspx. CUDA itself will be installed
    automatically when you install :code:`openmm` in the next steps.
  * If you have an AMD GPU and are using Linux or Windows, download the latest
    version of the drivers from https://support.amd.com.  To use the HIP
    platform (recommended), you also need to install HIP/ROCm by following the
    instructions at https://rocm.docs.amd.com.
  * On macOS, OpenCL is included with the operating system and is supported on
    macOS 10.10.3 or later.

3.a. To install with conda, open a command line terminal and type the following command
::

    conda install -c conda-forge openmm

With recent :code:`conda` versions (v4.8.4+), this will install a version of
OpenMM compiled with the latest version of CUDA supported by your drivers.
Alternatively you can request a version that is compiled for a specific CUDA
version with the command
::

    conda install -c conda-forge openmm cuda-version=12

where :code:`12` should be replaced with the particular CUDA version
you want to target.  We build packages for CUDA 11 and above.  Because different
CUDA releases are not binary compatible with each other, OpenMM can only work
with the particular CUDA version it was compiled with.

The conda package does not include the HIP platform.  If you have a recent AMD
GPU, we recommend installing with pip instead.

3.b. To install with pip, open a command line terminal and type the following command
::

    pip install openmm

The package installed with that command includes the OpenCL, CPU, and Reference
platforms.  To also install the CUDA platform (recommended if you have an NVIDIA
GPU), type
::

    pip install openmm[cuda12]

This will install a copy of the CUDA platform compiled with CUDA 12.  Alternatively,
if you have an AMD GPU, use this command to include the HIP platform (compiled
with HIP version 6)
::

    pip install openmm[hip6]

4. Verify your installation by typing the following command:
::

    python -m openmm.testInstallation

This command confirms that OpenMM is installed, checks whether GPU acceleration
is available (via the CUDA, OpenCL, and/or HIP platforms), and verifies that all
platforms produce consistent results.

