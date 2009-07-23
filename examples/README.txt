----------------------------------------------------------------
OpenMM(tm) Example programs and prototype C and Fortran wrappers
for OpenMM Preview Release 3, June, 2009.

See https://simtk.org/home/openmm.
For help go to the Advanced/Public Forums tab and post to the
"help" forum there.

OpenMM is part of SimTK, the Simbios simulation Toolkit, from 
the Simbios National Center for Physics-based Simulation of 
Biological Structures funded by the National Institutes of Health
through the NIH Roadmap for Medical Research, Grant U54 GM072970. 
Information on the National Centers for Biomedical Computing can 
be obtained here: http://nihroadmap.nih.gov/bioinformatics/.

OpenMM is developed under the supervision of Simbios P.I. Vijay 
Pande at Stanford. Any work that uses OpenMM should cite the 
following paper: M. S. Friedrichs, P. Eastman, V. Vaidyanathan, 
M. Houston, S. LeGrand, A. L. Beberg, D. L. Ensign, C. M. Bruns, 
V. S. Pande. “Accelerating Molecular Dynamic Simulation on 
Graphics Processing Units.” J. Comp. Chem., (2009), 30(6):864-872.

These simple "hello world" examples were developed by
Christopher Bruns and Michael Sherman for the OpenMM
workshop on June 24, 2009 at Stanford.
----------------------------------------------------------------


In addition to the example programs here, there are also two sets 
of wrappers to enable C and Fortran 95 to call the 
OpenMM API. However, even if your calling code is in C or Fortran, 
please consider instead writing your OpenMM-calling routines in 
C++, using 'extern "C"' so that your main program can call them. 
If that isn't feasible for you, try out our wrappers 
and let us know what troubles or better ideas you have. (Post to 
the OpenMM Forum at the above URL.)

Building the examples
---------------------

This directory includes a Makefile suitable for building the 
examples under Mac or Linux and possibly Cygwin, using the
gcc compiler suite. See MakefileNotes.txt for more info.

There is a file "NMakefile" which can be used with Microsoft's
NMake program to build the examples using Microsoft's "cl" compiler
for C++ and C, and Intel's "ifort" Fortran compiler.

There is a subdirectory here providing Visual Studio "solutions"
for building HelloArgon in C++, C, and Fortran. You will have to 
make your own for the other examples or just substitute a different 
source file for HelloArgon.cpp.

HelloArgon (C++, C, Fortran 95)
-------------------------------

This is the simplest example we could come up with that
does anything interesting. It is three argon atoms interacting
in a vacuum via van der Waals forces. It is primarily so you can
check that you are able to compile, link, and run correctly with
OpenMM. You will also get a many-frame pdb file generated
which you can view as an animation in VMD or most other
molecular viewers.

HelloSodiumChloride (C++, C, Fortran 95)
----------------------------------------

This example shows how we recommend using OpenMM so that you
can call it from an existing Molecular Dynamics code with
minimal disruption. The example contains Coulomb and van der 
Waals interactions and implicit solvation in a constant 
temperature simulation.

HelloEthane (C++ only)
----------------------

This example shows how to convey bond information to 
OpenMM. It is organized similarly to HelloSodiumChloride.

HelloWaterBox (C++ only)
------------------------

This example shows use of explicit solvent in a periodic box.
It is organized like the previous two.


C Wrapper
---------

This consists of two files:
    OpenMMCWrapper.h
    OpenMMCWrapper.cpp

Your C code (and the C examples above) include the header file
and link with OpenMMCWrapper.o along with the OpenMM library.
Consult the code or the example programs to figure out how
to use the wrappers; they are for the most part a very 
straightforward rehashing of the well-documented OpenMM C++
API.

Fortran Wrapper
---------------

The files that define the OpenMM Fortran bindings are:
     OpenMMFortranModule.f90
     OpenMMFortranWrapper.cpp

This consists of a Fortran 95 Module "OpenMM" so that your
Fortran program units have "use OpenMM" statements. You will
need access to the OpenMM*.mod files generated for the modules,
and you must also link with the OpenMMCWrapper.o as described
above for C, PLUS the OpenMMFortranWrapper.o. That's because 
the Fortran Wrapper makes use of the C Wrapper.









