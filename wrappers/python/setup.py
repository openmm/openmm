#

"""
setup.py: Used for building python wrappers for Simbios' OpenMM library.
"""
__author__ = "Randall J. Radmer"
__version__ = "1.0"


import os, sys, platform, glob, shutil
import struct

from distutils.core import setup

MAJOR_VERSION_NUM='4'
MINOR_VERSION_NUM='0'
BUILD_INFO='0'

def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)

def removeRecursive(dir):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isdir(path):
            removeRecursive(path)
        else:
            os.remove(path)
    os.rmdir(dir)

def removePackage(mod, verbose):
        try:
            pathList = mod.__path__
        except AttributeError:
            return
        if len(pathList) > 1:
           raise Exception("more than one item in simtk.__path__")
        simtkInstallPath = pathList[0]
        if os.path.exists(simtkInstallPath):
            if verbose:
                sys.stdout.write('REMOVING "%s"\n' % simtkInstallPath)
            removeRecursive(simtkInstallPath)

def uninstall(verbose=True):
    save_path=sys.path[:]
    sys.path=[]
    for item in save_path:
        if item!='.' and item!=os.getcwd():
            sys.path.append(item)
    try:
        import simtk.openmm as openmm
        removePackage(openmm, verbose)
    except ImportError:
        pass

    try:
        import simtk.unit as unit
        removePackage(unit, verbose)
    except ImportError:
        pass
    sys.path=save_path


def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):
    from distutils.core import Extension
    setupKeywords = {}
    setupKeywords["name"]              = "OpenMM"
    setupKeywords["version"]           = "%s.%s.%s" % (major_version_num,
                                                       minor_version_num,
                                                       build_info)
    setupKeywords["author"]            = "Randall J. Radmer"
    setupKeywords["author_email"]      = "radmer@stanford.edu"
    setupKeywords["license"]           = \
    "Python Software Foundation License (BSD-like)"
    setupKeywords["url"]               = "https://simtk.org/home/openmm"
    setupKeywords["download_url"]      = "https://simtk.org/home/openmm"
    setupKeywords["packages"]          = ["simtk",
                                          "simtk.unit",
                                          "simtk.openmm",
                                          "simtk.openmm.app",
                                          "simtk.openmm.app.internal"]
    setupKeywords["data_files"]        = []
    setupKeywords["package_data"]      = {"simtk" : [],
                                          "simtk.unit" : [],
                                          "simtk.openmm" : [],
                                          "simtk.openmm.app" : ['data/*.xml'],
                                          "simtk.openmm.app.internal" : []}
    setupKeywords["platforms"]         = ["Linux", "Mac OS X", "Windows"]
    setupKeywords["description"]       = \
    "Python wrapper for OpenMM (a C++ MD package)"
    setupKeywords["long_description"]  = \
    """OpenMM is a library which provides tools for modern molecular
    modeling simulation. As a library it can be hooked into any code,
    allowing that code to do molecular modeling with minimal extra
    coding (https://simtk.org/home/openmm).  This Python package
    gives access to the OpenMM API.
    """

    define_macros = [('MAJOR_VERSION', major_version_num),
                     ('MINOR_VERSION', minor_version_num)]

    libraries=['OpenMM',
               'OpenMMSerialization', 
               'OpenMMFreeEnergy',
               'OpenMMAmoeba',
              ]
    if 'OPENMM_USE_DEBUG_LIBS' in os.environ:
        if platform.system() == "Windows":
            raise Exception("use of OpenMM debug libs not supported on Win OS")
        else:
            sys.stdout.write("WARNING: using debug libs:\n")
            for ii in range(len(libraries)):
                libraries[ii]="%s_d" % libraries[ii]
                sys.stdout.write("%s\n" % libraries[ii])
            
    openmm_include_path = os.getenv('OPENMM_INCLUDE_PATH')
    if not openmm_include_path:
        reportError("Set OPENMM_INCLUDE_PATH to point to the include directory for OpenMM")
    openmm_lib_path = os.getenv('OPENMM_LIB_PATH')
    if not openmm_lib_path:
        reportError("Set OPENMM_LIB_PATH to point to the lib directory for OpenMM")

    extra_compile_args=[]
    extra_link_args=[]
    if platform.system() == "Windows":
        define_macros.append( ('WIN32', None) )
        define_macros.append( ('_WINDOWS', None) )
        define_macros.append( (' _MSC_VER', None) )
        extra_compile_args.append('/EHsc')
    else:
        if platform.system() == 'Darwin':
            macVersion = [int(x) for x in platform.mac_ver()[0].split('.')]
            if tuple(macVersion) < (10, 6):
                os.environ['MACOSX_DEPLOYMENT_TARGET']='10.5'
            extra_link_args.append('-Wl,-rpath,@loader_path/OpenMM')


    library_dirs=[openmm_lib_path]
    include_dirs=openmm_include_path.split(';')

    setupKeywords["ext_modules"] = [
       Extension(name = "simtk.openmm._openmm",
                 sources = ["src/swig_doxygen/OpenMMSwig.cxx"],
                 include_dirs = include_dirs,
                 define_macros = define_macros,
                 library_dirs = library_dirs,
                 libraries = libraries,
                 extra_compile_args=extra_compile_args,
                 extra_link_args=extra_link_args,
                 )
    ]

    outputString = ''
    firstTab     = 40
    secondTab    = 60
    for key in sorted( setupKeywords.iterkeys() ):
         value         = setupKeywords[key]
         outputString += key.rjust(firstTab) + str( value ).rjust(secondTab) + "\n"
    
    print "%s" % outputString

    return setupKeywords
    

def main():
    if sys.version_info < (2, 6):
        reportError("OpenMM requires Python 2.6 or better.")
    if sys.version_info >= (3,):
        reportError("OpenMM has not been tested with Python 3.0 or higher.")
    if platform.system() == 'Darwin':
        macVersion = [int(x) for x in platform.mac_ver()[0].split('.')]
        if tuple(macVersion) < (10, 5):
            reportError("OpenMM requires Mac OS X Leopard (10.5) or better.")
    uninstall()
    setupKeywords=buildKeywordDictionary()
    setup(**setupKeywords)

if __name__ == '__main__':
    main()


