This swig based, wrapper code only works from a Bash shell.

You *must* have the following installed:
1) Python 2.5 or 2.6
2) swig 2.0.0 or better (earlier versions are likely to fail)
3) Some type of XSLT processor (xsltproc, saxon, etc)
4) Doxygen (tested with version 1.4.7, but others should also work)
5) py-dom-xpath (tested with version 0.1, but others should also work):
   http://code.google.com/p/py-dom-xpath


To build the wrapper codes, the following must be done
1) Run doxygen to generate xml files containing OpenMM API info (doxygen
   uses the OpenMM header files to do this).
2) Run a python script that builds SWIG input files based on info in the
   doxygen xml files
3) Run SWIG to generate python and C code that wraps the OpenMM libraries

We have created a single Bash script that does the full build for you.
But before running it, you must create and export environment variables
containing the OpenMM lib path (OPENMM_LIB_PATH) and the include directory
path (OPENMM_INCLUDE_PATH).  The include directory should contain the "OpenMM.h"
file and an "openmm" directory.  The lib directory should contain
the OpenMM dynamic library (OpenMM.dll, libOpenMM.so or libOpenMM.dylib) and
plugins.  For example, when I build OpenMM on a Bash system, I set the
install dir to $HOME (in ccmake). Then I type the following two lines at
the command line:
export OPENMM_LIB_PATH=$HOME/lib
export OPENMM_INCLUDE_PATH=$HOME/include

Next run the Bash build script as follows:
/bin/bash buildSwigWrapper.sh

You may need to edit this script to work on your system

If all went well, the wrapper code is now updated; cd ../.. and follow
instuctions for building the PyOpenMM libraries.


1. Edit swigInputConfig.py to add new classes, ...
2. Add template <class T> to XmlSerializer (serialize/deserialize)

    class XmlSerializer {
    public:
    
       %apply std::ostream & OUTPUT { std::ostream & stream };
       template <class T>static void serialize(const T *object, const std::string
    &rootName, std::ostream &stream) ;
       %clear std::ostream & stream;
       %apply std::istream & OUTPUT { std::istream & stream };
       template <class T>static T* deserialize(std::istream &stream) ;
       %clear std::istream & stream;
    };



