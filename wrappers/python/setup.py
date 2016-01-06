from __future__ import print_function
import re
import os
import sys
import shutil
import platform

if 'bdist_wheel' in sys.argv:
    from setuptools import setup, Extension
else:
    from distutils.core import setup, Extension

from distutils.command.install_lib import install_lib as _install_lib

# Globals configured by CMake
MAJOR_VERSION_NUM = '@OPENMM_MAJOR_VERSION@'
MINOR_VERSION_NUM = '@OPENMM_MINOR_VERSION@'
BUILD_INFO = '@OPENMM_BUILD_VERSION@'
CMAKE_SOURCE_DIR = '@CMAKE_SOURCE_DIR@'
ALL_PLUGIN_TARGET_LOCATIONS = list(filter(None, '@ALL_PLUGIN_TARGET_LOCATIONS@'.split(';')))
ALL_LIBRARIES = list(filter(None, '@ALL_LIBRARIES@'.split(';')))
ALL_LIBRARY_LOCATIONS = list(filter(None, '@ALL_LIBRARY_LOCATIONS@'.split(';')))
WRAPPER_INCLUDE_DIRS =  list(filter(None, '@WRAPPER_INCLUDE_DIRS@'.split(';')))

__author__ = "Peter Eastman"
__version__ = "%s.%s" % (MAJOR_VERSION_NUM, MINOR_VERSION_NUM)


def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)


class install_lib(_install_lib):
    def install(self):
        self.copyOpenMMLibraries()
        self.copyOpenMMPlugins()
        self.copyLicenses()
        return _install_lib.install(self)

    def copyOpenMMLibraries(self):
        dest = os.path.join(self.build_dir, 'simtk/openmm/')
        if not os.path.exists(dest):
            os.makedirs(dest)
        for src in ALL_LIBRARY_LOCATIONS:
            shutil.copy2(src, dest)

    def copyOpenMMPlugins(self):
        dest = os.path.join(self.build_dir, 'simtk/openmm/plugins/')
        if os.path.exists(dest):
            shutil.rmtree(dest)
        if not os.path.exists(dest):
            os.makedirs(dest)

        for src in ALL_PLUGIN_TARGET_LOCATIONS:
            shutil.copy2(src, dest)

    def copyLicenses(self):
        dest = os.path.join(self.build_dir, 'simtk/openmm/licenses')
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(os.path.join(CMAKE_SOURCE_DIR, 'docs-source/licenses'), dest)


def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):

    setupKeywords = {
        "name": "OpenMM",
        "version":
        "%s.%s.%s" % (major_version_num, minor_version_num, build_info),
        "author": "Peter Eastman",
        "license": "LGPL",
        "url": "https://simtk.org/home/openmm",
        "download_url": "https://simtk.org/home/openmm",
        "packages": ["simtk", "simtk.unit", "simtk.openmm", "simtk.openmm.app",
                     "simtk.openmm.app.internal",
                     "simtk.openmm.app.internal.charmm",
                     "simtk.openmm.app.internal.pdbx",
                     "simtk.openmm.app.internal.pdbx.reader",
                     "simtk.openmm.app.internal.pdbx.writer"],
        "data_files": [],
        "package_data": {"simtk": [],
                         "simtk.unit": [],
                         "simtk.openmm": [],
                         "simtk.openmm.app": ['data/*.xml', 'data/*.pdb'],
                         "simtk.openmm.app.internal": []},
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
    if 'setuptools' in sys.modules:
        setupKeywords['zip_safe'] = False

    define_macros = [('MAJOR_VERSION', major_version_num),
                     ('MINOR_VERSION', minor_version_num)]

    extra_compile_args = []
    extra_link_args = []
    if platform.system() == "Windows":
        define_macros.append(('WIN32', None))
        define_macros.append(('_WINDOWS', None))
        define_macros.append((' _MSC_VER', None))
        extra_compile_args.append('/EHsc')
    elif platform.system() == 'Darwin':
        extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7',
                            '-Wl,-rpath,@loader_path/.']
        # Hard-code CC and CXX to clang, since gcc/g++ will *not* work with
        # Anaconda, despite the fact that distutils will try to use them.
        # System Python, homebrew, and MacPorts on Macs will always use
        # clang, so this hack should always work and fix issues with users
        # that have GCC installed from MacPorts or homebrew *and* Anaconda
        os.environ['CC'] = 'clang'
        os.environ['CXX'] = 'clang++'
    elif platform.system() == 'Linux':
        extra_link_args = ['-Wl,-rpath,$ORIGIN/.', ]

    library_dirs = list(set([os.path.dirname(lib) for lib in ALL_LIBRARY_LOCATIONS]))

    extensionArgs = {"name": "simtk.openmm._openmm",
                     "sources": ["src/swig_doxygen/OpenMMSwig.cxx"],
                     "include_dirs": WRAPPER_INCLUDE_DIRS,
                     "define_macros": define_macros,
                     "libraries": ALL_LIBRARIES,
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
    setup(**setupKeywords)


if __name__ == '__main__':
    main()
