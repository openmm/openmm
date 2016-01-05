"""
setup.py: Used for building python wrappers for Simbios' OpenMM library.
"""
import ast
import re
import os
import sys
import shutil
import platform

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from distutils.command.install_lib import install_lib as _install_lib

MAJOR_VERSION_NUM = '@OPENMM_MAJOR_VERSION@'
MINOR_VERSION_NUM = '@OPENMM_MINOR_VERSION@'
BUILD_INFO = '@OPENMM_BUILD_VERSION@'
IS_RELEASED = False

__author__ = "Peter Eastman"
__version__ = "%s.%s" % (MAJOR_VERSION_NUM, MINOR_VERSION_NUM)


def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)


def libraryFilename(libName):
    if platform.system() == 'Windows':
        return '%s.dll' % libName
    elif platform.system() == 'Darwin':
        return 'lib%s.dylib' % libName
    elif platform.system() == 'Linux':
        return 'lib%s.so' % libName
    else:
        raise ValueError('Unknown platform: "%s"' % platform.system())


def getLibraries():
    libraries = ['OpenMM',
                 'OpenMMAmoeba',
                 'OpenMMRPMD',
                 'OpenMMDrude']

    if 'OPENMM_USE_DEBUG_LIBS' in os.environ:
        if platform.system() == "Windows":
            raise Exception("use of OpenMM debug libs not supported on Win OS")
        else:
            sys.stdout.write("WARNING: using debug libs:\n")
            for ii in range(len(libraries)):
                libraries[ii] = "%s_d" % libraries[ii]
                sys.stdout.write("%s\n" % libraries[ii])
    return libraries


class install_lib(_install_lib):
    def install(self):
        openMMLibPath = os.getenv('OPENMM_LIB_PATH')
        if not openMMLibPath:
            reportError("Set OPENMM_LIB_PATH to point to the lib directory for OpenMM")

        self.copyOpenMMLibraries(openMMLibPath)
        self.copyOpenMMPlugins(openMMLibPath)
        return _install_lib.install(self)

    def copyOpenMMLibraries(self, openMMLibPath):
        dest = os.path.join(self.build_dir, 'simtk/openmm')
        if not os.path.exists(dest):
            os.makedirs(dest)
        for libName in getLibraries():
            libFileName = libraryFilename(libName)
            src = os.path.join(openMMLibPath, libFileName)
            shutil.copy2(src, dest)

    def copyOpenMMPlugins(self, openMMLibPath):
        dest = os.path.join(self.build_dir, 'simtk/openmm/plugins')
        if os.path.exists(dest):
            shutil.rmtree(dest)

        src = os.path.join(openMMLibPath, 'plugins')
        shutil.copytree(src, dest)


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
                       'isrelease': str(IS_RELEASED)})
    finally:
        a.close()


def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):


    setupKeywords = {
        "name": "OpenMM",
        "version": "%s.%s.%s" % (major_version_num, minor_version_num, build_info),
        "author":  "Peter Eastman",
        "license": "LGPL",
        "url": "https://simtk.org/home/openmm",
        "download_url": "https://simtk.org/home/openmm",
        "packages": ["simtk",
                     "simtk.unit",
                     "simtk.openmm",
                     "simtk.openmm.app",
                     "simtk.openmm.app.internal",
                     "simtk.openmm.app.internal.charmm",
                     "simtk.openmm.app.internal.pdbx",
                     "simtk.openmm.app.internal.pdbx.reader",
                     "simtk.openmm.app.internal.pdbx.writer"],
        "data_files": [],
        "package_data": {"simtk" : [],
                         "simtk.unit" : [],
                         "simtk.openmm" : [],
                         "simtk.openmm.app" : ['data/*.xml', 'data/*.pdb'],
                         "simtk.openmm.app.internal" : []},
        "platforms": ["Linux", "Mac OS X", "Windows"],
        "description": "High performance molecular simulation on GPUs",
        "long_description":
            """OpenMM is a toolkit for molecular simulation. It can be used either as a
            stand-alone application for running simulations, or as a library you call
            from your own code. It provides a combination of extreme flexibility
            (through custom forces and integrators), openness, and high performance
            (especially on recent GPUs) that make it truly unique among simulation codes.
            """,
        "cmdclass": {"install_lib": install_lib},
    }

    define_macros = [('MAJOR_VERSION', major_version_num),
                     ('MINOR_VERSION', minor_version_num)]


    openmm_include_path = os.getenv('OPENMM_INCLUDE_PATH')
    if not openmm_include_path:
        reportError("Set OPENMM_INCLUDE_PATH to point to the include directory for OpenMM")

    extra_compile_args=[]
    extra_link_args=[]
    if platform.system() == "Windows":
        define_macros.append( ('WIN32', None) )
        define_macros.append( ('_WINDOWS', None) )
        define_macros.append( (' _MSC_VER', None) )
        extra_compile_args.append('/EHsc')
    elif platform.system() == 'Darwin':
        extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl,-rpath,@loader_path/.']
        # Hard-code CC and CXX to clang, since gcc/g++ will *not* work with
        # Anaconda, despite the fact that distutils will try to use them.
        # System Python, homebrew, and MacPorts on Macs will always use
        # clang, so this hack should always work and fix issues with users
        # that have GCC installed from MacPorts or homebrew *and* Anaconda
        os.environ['CC'] = 'clang'
        os.environ['CXX'] = 'clang++'
    elif platform.system() == 'Linux':
        extra_link_args = ['-Wl,-rpath,$ORIGIN/.',]

    libraries = getLibraries()
    library_dirs = [os.getenv('OPENMM_LIB_PATH')]
    include_dirs = openmm_include_path.split(';')

    extensionArgs = {"name": "simtk.openmm._openmm",
                     "sources": ["src/swig_doxygen/OpenMMSwig.cxx"],
                     "include_dirs": include_dirs,
                     "define_macros": define_macros,
                     "libraries": libraries,
                     "library_dirs": library_dirs,
                     "extra_compile_args": extra_compile_args,
                     "extra_link_args": extra_link_args}

    setupKeywords["ext_modules"] = [Extension(**extensionArgs)]
    return setupKeywords


def main():
    if sys.version_info < (2, 6):
        reportError("OpenMM requires Python 2.6 or better.")
    if platform.system() == 'Darwin':
        macVersion = [int(x) for x in platform.mac_ver()[0].split('.')]
        if tuple(macVersion) < (10, 5):
            reportError("OpenMM requires Mac OS X Leopard (10.5) or better.")

    setupKeywords = buildKeywordDictionary()
    writeVersionPy()

    setup(**setupKeywords)


if __name__ == '__main__':
    main()


