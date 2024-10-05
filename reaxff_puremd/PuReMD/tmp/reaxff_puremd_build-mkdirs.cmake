# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD")
  file(MAKE_DIRECTORY "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD")
endif()
file(MAKE_DIRECTORY
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src/reaxff_puremd_build-build"
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD"
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/tmp"
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src/reaxff_puremd_build-stamp"
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src"
  "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src/reaxff_puremd_build-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src/reaxff_puremd_build-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/babaid/repos/openmm-puremd/reaxff_puremd/PuReMD/src/reaxff_puremd_build-stamp${cfgdir}") # cfgdir has leading slash
endif()
