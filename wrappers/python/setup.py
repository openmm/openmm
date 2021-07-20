"""
setup.py: Used for building python wrappers for Simbios' OpenMM library.
"""
import ast
import re
import os
import sys
import platform
import numpy
from setuptools import setup
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
           raise Exception("more than one item in openmm.__path__")
        installPath = pathList[0]
        if os.path.exists(installPath):
            if verbose:
                sys.stdout.write('REMOVING "%s"\n' % installPath)
            removeRecursive(installPath)

def uninstall(verbose=True):
    save_path=sys.path[:]
    sys.path=[]
    for item in save_path:
        if item!='.' and item!=os.getcwd():
            sys.path.append(item)
    try:
        import simtk.openmm
        removePackage(simtk.openmm, verbose)
    except ImportError:
        pass

    try:
        import simtk.unit as unit
        removePackage(unit, verbose)
    except ImportError:
        pass

    try:
        import openmm
        removePackage(openmm, verbose)
    except ImportError:
        pass
    sys.path=save_path


def writeVersionPy(filename="openmm/version.py", major_version_num=MAJOR_VERSION_NUM,
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
%(path)s

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

    if 'OPENMM_MAKE_WHEEL' in os.environ:
        path = '''
import os as _os
_base_path = _os.path.abspath(_os.path.dirname(__file__))
openmm_library_path = _os.path.join(_base_path, 'lib')
openmm_include_path = _os.path.join(_base_path, 'include')
        '''
    else:
        path = '''
openmm_library_path = '{}'
openmm_include_path = '{}'
'''.format(os.getenv('OPENMM_LIB_PATH'),os.getenv('OPENMM_INCLUDE_PATH'))

    with open(filename, 'w') as a:
        a.write(cnt % {'version': version,
                       'full_version' : full_version,
                       'git_revision' : git_revision,
                       'isrelease': str(IS_RELEASED),
                       'path': path})

def lib_and_include_files(base_dir):
    plat = platform.system()
    if plat == 'Darwin':
        lib_ext = '*.dylib'
    elif plat == 'Linux':
        lib_ext = "*.so"
    else:
        lib_ext = "*.dll"
    lib_dir = os.path.join(base_dir, 'lib')
    inc_dir = os.path.join(base_dir, 'include')
    paths = []
    for (path, directories, filenames) in os.walk(lib_dir):
        unique_bases = set(f.split('.')[0] for f in filenames)
        for file_base in unique_bases:
            matching = [f for f in filenames if f.split('.')[0] == file_base]
            if len(matching)==1:
                keep = matching[0]
            else:
                # Find the most fully versioned file and keep that
                matching = list(sorted(matching, key=lambda f:len(f.split('.'))))
                keep = matching[-1]
            paths.append(os.path.join(os.path.relpath(path, base_dir), keep))
    for (path, directories, _) in os.walk(inc_dir):
        paths.append(os.path.join(os.path.relpath(path, base_dir), '*.h'))
    return paths

def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):
    from distutils.core import Extension
    setupKeywords = {}
    setupKeywords["name"]              = "OpenMM"
    setupKeywords["version"]           = "%s.%s.%s" % (major_version_num,
                                                       minor_version_num,
                                                       build_info)
    setupKeywords["author"]            = "Peter Eastman"
    setupKeywords["license"]           = \
    "Python Software Foundation License (BSD-like)"
    setupKeywords["url"]               = "https://openmm.org"
    setupKeywords["download_url"]      = "https://openmm.org"
    setupKeywords["packages"]          = [
                                          "simtk",
                                          "simtk.unit",
                                          "simtk.openmm",
                                          "simtk.openmm.app",
                                          "openmm",
                                          "openmm.unit",
                                          "openmm",
                                          "openmm.app",
                                          "openmm.app.internal",
                                          "openmm.app.internal.charmm",
                                          "openmm.app.internal.pdbx",
                                          "openmm.app.internal.pdbx.reader",
                                          "openmm.app.internal.pdbx.writer"]
    setupKeywords["data_files"]        = []
    setupKeywords["package_data"]      = {"openmm" : [],
                                          "openmm.app" : ['data/*.xml', 'data/*.pdb', 'data/amber14/*.xml', 'data/charmm36/*.xml'],
                                          "openmm.app.internal" : []}
    if 'OPENMM_MAKE_WHEEL' in os.environ:
        setupKeywords['package_data']['openmm'].extend(lib_and_include_files(os.path.join(os.curdir,'openmm')))
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
    

    extra_compile_args=['-std=c++11']
    extra_link_args=[]
    if platform.system() == "Windows":
        define_macros.append( ('WIN32', None) )
        define_macros.append( ('_WINDOWS', None) )
        define_macros.append( (' _MSC_VER', None) )
        extra_compile_args.append('/EHsc')
    else:
        if platform.system() == 'Darwin':
            extra_compile_args += ['-stdlib=libc++']
            if 'OPENMM_MAKE_WHEEL' in os.environ:
                rpath = os.path.join('@loader_path','lib')
            else:
                rpath = openmm_lib_path
            extra_link_args += ['-stdlib=libc++', '-Wl', '-rpath', rpath]
            if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ and platform.processor() != 'arm':
                extra_compile_args += ['-mmacosx-version-min=10.7']
                extra_link_args += ['-mmacosx-version-min=10.7']
            if 'OPENMM_MAKE_WHEEL' not in os.environ:
                # Hard-code CC and CXX to clang, since gcc/g++ will *not* work with
                # Anaconda, despite the fact that distutils will try to use them.
                # System Python, homebrew, and MacPorts on Macs will always use
                # clang, so this hack should always work and fix issues with users
                # that have GCC installed from MacPorts or homebrew *and* Anaconda
                if 'CC' not in os.environ:
                    os.environ['CC'] = 'clang'
                    os.environ['CXX'] = 'clang++'
        elif 'OPENMM_MAKE_WHEEL' in os.environ:
            extra_link_args += ['-Wl,-rpath,'+os.path.join('$ORIGIN','lib')]

    library_dirs=[openmm_lib_path]
    include_dirs=openmm_include_path.split(';')
    include_dirs.append(numpy.get_include())

    extensionArgs = {"name": "openmm._openmm",
                    "sources": ["src/swig_doxygen/OpenMMSwig.cxx"],
                    "include_dirs": include_dirs,
                    "define_macros": define_macros,
                    "library_dirs": library_dirs,
                    "libraries": libraries,
                    "extra_compile_args": extra_compile_args,
                    "extra_link_args": extra_link_args}
    if platform.system() != "Windows" and 'OPENMM_MAKE_WHEEL' not in os.environ:
        extensionArgs["runtime_library_dirs"] = library_dirs
    setupKeywords["ext_modules"] = [Extension(**extensionArgs)]
    setupKeywords["ext_modules"] += cythonize('openmm/app/internal/*.pyx')

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
    import os
    if 'OPENMM_MAKE_WHEEL' not in os.environ:
        try:
            uninstall()
        except:
            pass
    setupKeywords=buildKeywordDictionary()
    writeVersionPy()
    setup(**setupKeywords)

if __name__ == '__main__':
    main()
