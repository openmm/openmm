"""
setup.py: Used for building python wrappers for Simbios' OpenMM library.
"""
import ast
import re
import os
import sys
import platform
import numpy
from distutils.core import setup
from Cython.Build import cythonize

MAJOR_VERSION_NUM='@OPENMM_MAJOR_VERSION@'
MINOR_VERSION_NUM='@OPENMM_MINOR_VERSION@'
BUILD_INFO='@OPENMM_BUILD_VERSION@'
IS_RELEASED = False

__author__ = "Peter Eastman"
__version__ = "%s.%s" % (MAJOR_VERSION_NUM, MINOR_VERSION_NUM)

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


def writeVersionPy(filename="simtk/openmm/version.py", major_version_num=MAJOR_VERSION_NUM,
                     minor_version_num=MINOR_VERSION_NUM, build_info=BUILD_INFO):
    """Write a version.py file into the python source directory before installation.
    If a version.py file already exists, we assume that it contains only the git_revision
    information, since from within this python session in the python staging directory, we're
    not in the version controlled directory hierarchy.

    When cmake is copying files into the PYTHON_STAGING_DIRECTORY, it will write the
    git revision to version.py. We read that, and then overwrite it.
    """

    cnt = """
# THIS FILE IS GENERATED FROM OPENMM SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s
openmm_library_path = r'%(path)s'

if not release:
    version = full_version
"""

    if os.path.exists(filename):
        # git_revision is written to the file by cmake
        with open(filename) as f:
            text = f.read()
            match = re.search(r"git_revision\s+=\s+(.*)", text, re.MULTILINE)
        try:
            git_revision = ast.literal_eval(match.group(1))
        except:
            # except anything, including no re match or
            # literal_eval failing
            git_revision = 'Unknown'
    else:
        git_revision = 'Unknown'

    version = full_version = '%s.%s.%s' % (major_version_num, minor_version_num, build_info)
    if not IS_RELEASED:
        full_version += '.dev-' + git_revision[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': version,
                       'full_version' : full_version,
                       'git_revision' : git_revision,
                       'isrelease': str(IS_RELEASED),
                       'path': os.getenv('OPENMM_LIB_PATH')})
    finally:
        a.close()


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
                                          "simtk.openmm.app.internal",
                                          "simtk.openmm.app.internal.charmm",
                                          "simtk.openmm.app.internal.pdbx",
                                          "simtk.openmm.app.internal.pdbx.reader",
                                          "simtk.openmm.app.internal.pdbx.writer"]
    setupKeywords["data_files"]        = []
    setupKeywords["package_data"]      = {"simtk" : [],
                                          "simtk.unit" : [],
                                          "simtk.openmm" : [],
                                          "simtk.openmm.app" : ['data/*.xml', 'data/*.pdb', 'data/amber14/*.xml', 'data/charmm36/*.xml'],
                                          "simtk.openmm.app.internal" : []}
    setupKeywords["platforms"]         = ["Linux", "Mac OS X", "Windows"]
    setupKeywords["description"]       = \
    "Python wrapper for OpenMM (a C++ MD package)"
    setupKeywords["long_description"]  = \
    """OpenMM is a toolkit for molecular simulation. It can be used either as a
    stand-alone application for running simulations, or as a library you call
    from your own code. It provides a combination of extreme flexibility
    (through custom forces and integrators), openness, and high performance
    (especially on recent GPUs) that make it truly unique among simulation codes.
    """

    define_macros = [('MAJOR_VERSION', major_version_num),
                     ('MINOR_VERSION', minor_version_num)]

    libraries=['OpenMM',
               'OpenMMAmoeba',
               'OpenMMRPMD',
               'OpenMMDrude',
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
            extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
            extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_lib_path]
            # Hard-code CC and CXX to clang, since gcc/g++ will *not* work with
            # Anaconda, despite the fact that distutils will try to use them.
            # System Python, homebrew, and MacPorts on Macs will always use
            # clang, so this hack should always work and fix issues with users
            # that have GCC installed from MacPorts or homebrew *and* Anaconda
            os.environ['CC'] = 'clang'
            os.environ['CXX'] = 'clang++'

    library_dirs=[openmm_lib_path]
    include_dirs=openmm_include_path.split(';')
    include_dirs.append(numpy.get_include())

    extensionArgs = {"name": "simtk.openmm._openmm",
                    "sources": ["src/swig_doxygen/OpenMMSwig.cxx"],
                    "include_dirs": include_dirs,
                    "define_macros": define_macros,
                    "library_dirs": library_dirs,
                    "libraries": libraries,
                    "extra_compile_args": extra_compile_args,
                    "extra_link_args": extra_link_args}
    if platform.system() != "Windows":
        extensionArgs["runtime_library_dirs"] = library_dirs
    setupKeywords["ext_modules"] = [Extension(**extensionArgs)]
    setupKeywords["ext_modules"] += cythonize('simtk/openmm/app/internal/*.pyx', language='c++')

    outputString = ''
    firstTab     = 40
    secondTab    = 60
    for key in sorted(iter(setupKeywords)):
         value         = setupKeywords[key]
         outputString += key.rjust(firstTab) + str( value ).rjust(secondTab) + "\n"

    sys.stdout.write("%s" % outputString)

    return setupKeywords


def main():
    if sys.version_info < (2, 7):
        reportError("OpenMM requires Python 2.7 or better.")
    if platform.system() == 'Darwin':
        macVersion = [int(x) for x in platform.mac_ver()[0].split('.')]
        if tuple(macVersion) < (10, 5):
            reportError("OpenMM requires Mac OS X Leopard (10.5) or better.")
    try:
        uninstall()
    except:
        pass
    setupKeywords=buildKeywordDictionary()
    writeVersionPy()
    setup(**setupKeywords)

if __name__ == '__main__':
    main()


