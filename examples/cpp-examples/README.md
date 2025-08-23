## OpenMM C++, C, and Fortran API Examples

These simple "hello world" examples were developed by Christopher Bruns and
Michael Sherman.  They demonstrate use of the C++, C, and Fortran APIs for
OpenMM.

### Example Details

Consult the OpenMM user guide for
[more information on these examples](https://docs.openmm.org/latest/userguide/library/03_tutorials.html),
as well as
[use of the C and Fortran API wrappers](https://docs.openmm.org/latest/userguide/library/05_languages_not_cpp.html).

#### `HelloArgon` (C++, C, Fortran 95)

This is the simplest example we could come up with that
does anything interesting. It is three argon atoms interacting
in a vacuum via van der Waals forces. It is primarily so you can
check that you are able to compile, link, and run correctly with
OpenMM. You will also get a many-frame PDB file generated
which you can view as an animation in VMD or most other
molecular viewers.

#### `HelloSodiumChloride` (C++, C, Fortran 95)

This example shows how we recommend using OpenMM so that you
can call it from an existing molecular dynamics code with
minimal disruption. The example contains Coulomb and van der
Waals interactions and implicit solvation in a constant
temperature simulation.

#### `HelloEthane` (C++ only)

This example shows how to convey bond information to
OpenMM. It is organized similarly to `HelloSodiumChloride`.

#### `HelloWaterBox` (C++ only)

This example shows use of explicit solvent in a periodic box.
It is organized like the previous two.

### Building the examples

This directory includes a Makefile suitable for building the  examples under Mac
or Linux and possibly Cygwin, using the GCC compiler suite.  There is also an
NMakefile for use with Microsoft's `nmake` program. For Fortran, the Makefile
expects `gfortran` to be in the path while NMakefile expects `ifort` (Intel
Fortran).  In the `VisualStudio` subdirectory, there are Visual Studio
"solutions" for building `HelloArgon` in C++, C, and Fortran. You will have to
make your own for the other examples or just substitute a different source file
for HelloArgon.

You must already have the OpenMM binaries installed to build the examples from
source.  See the
[user guide](https://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm)
for more information.  You may need to slightly edit the Makefile or NMakefile
to make  it run on your system, depending where you installed OpenMM and  the
particular requirements of your compiler versions for mixed Fortran/C++
programming.

Type `make` (or `make default`) to get just one C++ example built
(`HelloArgon`). Make sure it runs.  To compile all example programs type
`make all`. That includes Fortran examples though so you will see failures
unless you have `gfortran` installed (Linux and Mac) or `ifort` installed
(Windows). However, the C++ and C examples should compile with just `g++` or
`cl`.  To build just one example, type `make HelloArgonInC` or whatever.

Before you run the executables, remember to add the OpenMM dynamic library
directory to your library path.  The simplest way to do this is to type the
following commands:

- For linux (for the bash shell, assuming installation was done in the default
  location `/usr/local/openmm`):
  ```bash
  $ export LD_LIBRARY_PATH=/usr/local/openmm/lib
  ```
- For Mac (for the `bash` shell, assuming installation was done in the default
  location `/usr/local/openmm`):
  ```bash
  $ export DYLD_LIBRARY_PATH=/usr/local/openmm/lib
  ```
- For Windows (command tool, assuming installation was done in the default
  location `C:\Program Files\OpenMM`):
  ```bat
  C:\> set path=%path%;C:\Program Files\OpenMM\lib
  ```
