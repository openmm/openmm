===========================================
 PyOpenMM: Python wrappers for OpenMM
===========================================

URL: https://simtk.org/home/pyopenmm

PyOpenMM is a python API that allows access to OpenMM, which is a library
that provides tools for performing GPU accelerated molecular simulations.
See the OpenMM project home page for more details:
https://simtk.org/home/openmm

Note that you do *not* need to download the OpenMM package.  This
package includes all requried OpenMM libraries for Windows, Mac OS X,
Linux, and 64 bit Linux. However, if you wish to make use of the
GPU acceleration, you will need to have a supported GPU, and you will
need to install the associated CUDA Toolkit and drivers.


INTENDED AUDIENCE:
This README is intended for people who downloaded the PyOpenMM package
from the simtk.org website.


ADVANCED USERS:
if you are building or installing your own OpenMM libraries, you will also
need to rebuild the python wrappers.  Note that most users will not
need (or want) to do this.  If you really do need to rebuild the python
wrappers, change to the following subdirectory and read its README.txt file,
then return here to finish the installation:
cd src/swig_doxygen


REQUIREMENTS
1) Python 2.6 (Pythons 2.4 and 2.5 no longer supported; Python 2.7 might work)
2) One of the OS types supported by OpenMM:
  a) Windows XP or better
  b) Mac Leopard
  c) Snow Leopard, assuming 32-bit python. See the README_MAC_64bit.txt file.
  d) Linux (we test using CentOS 5, but others also should work)
  e) Linux64 (we test using CentOS 5, but others should work)
3) A C/C++ compiler.
  a) For Windows, we suggest Visual Studio 2008
  b) For Mac Leopard and Snow Leopard, you may need to install Apple's Xcode
     development environment (this will install otool and install_name_tool
     which might also be needed).
  c) For Linux, we use gcc-c++ 4.1


INSTALL
First unzip the package, change to the new directory, and run
the following two commands from a command line:
 python setup.py build
 python setup.py install


TEST INSTALL
Next run the following test to see if the install worked (it
should report "OK"):
 python test/test_openmm.py
This should return "OK"


USING A GPU
If you would like to make use of OpenMM's GPU acceleration, you
will have to have one of the supported GPU's, and you will need
to install the Cuda libraries.  To download the Cuda libraries,
check the NVIDIA website (try
http://www.nvidia.com/object/cuda_get.html).
For more help, see the OpenMM discussion forums at the following URL:
https://simtk.org/forum/?group_id=161)

Once the Cuda librareis are installed, python may need help finding
them.  One easy way to do this on Unix type machines is to set the
LD_LIBRARY_PATH to point to the Cuda installation.
On my Bash system, I could do the following at the command line:
 export LD_LIBRARY_PATH=/usr/local/cuda/lib
On a windows machine, adding the Cuda installation path to
your system PATH variable should work.

After installing you may need to reboot.

TEST GPU
If you have a GPU, run the following to see if PyOpenMM is accessing it:
 python scripts/pluginLoadingCheck.py
If you do not have a GPU, or it is not found, you should get the following:
 OpenMM found no plugins
If you have a GPU and it's found, you will see a list of GPU libraries that
have been loaded.  For example, on my system I see the following:
 OpenMM loaded the following plugin(s)
 $PYTHONPATH/simtk/chem/OpenMM/plugins/libOpenMMCuda.so

If you are unable to load the GPU (Cuda) plugins, check that you installed
the Cuda drivers and Toolkit correctly.


PYOPENMM USAGE EXAMPLES
Once you have tested the installation, switch to the example
directory (cd examples), and try running one of the simulations.
For example, run the argon simulation by changing to the examples/argon
directory and typing the following:
 python runArgon.py argon.pdb

If all goes well, this will produce a PDB formatted file called argon.pdb,
which can be visualized using a program such as VMD.


FEEDBACK
Please address all comments, criticisms, suggestions, etc to one of us:
Mark Friedrichs at friedrim@stanford.edu
Peter Eastman at cmbruns@stanford.edu 
Randy Radmer at radmer@stanford.edu 


