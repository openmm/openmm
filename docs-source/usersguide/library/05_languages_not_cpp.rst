.. _using-openmm-with-software-written-in-languages-other-than-c++:

Using OpenMM with Software Written in Languages Other than C++
##############################################################

Although the native OpenMM API is object-oriented C++ code, it is possible to
directly translate the interface so that it is callable from C, Fortran 95, and
Python with no substantial conceptual changes. We have developed a
straightforward mapping for these languages that, while perhaps not the most
elegant possible, has several advantages:

* Almost all documentation, training, forum discussions, and so on are equally
  useful to users of all these languages. There are syntactic differences of
  course, but all the important concepts remain unchanged.
* We are able to generate the C, Fortran, and Python APIs from the C++ API.
  Obviously, this reduces development effort, but more importantly it means that
  the APIs are likely to be error-free and are always available immediately when
  the native API is updated.
* Because OpenMM performs expensive operations “in bulk” there is no noticeable
  overhead in accessing these operations through the C, Fortran, or Python APIs.
* All symbols introduced to a C or Fortran program begin with the prefix
  “\ :code:`OpenMM_`\ ” so will not interfere with symbols already in use.


*Availability of APIs in other languages:*  All necessary C and Fortran
bindings are built in to the main OpenMM library; no separate library is
required.  The Python wrappers are contained in a module that is distributed
with OpenMM and that can be installed by executing its setup.py script in the
standard way.

(This doesn’t apply to most users: if you are building your own OpenMM from
source using CMake and want the API bindings generated, be sure to enable the
:code:`OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS` option for C and Fortran, or
:code:`OPENMM_BUILD_PYTHON_WRAPPERS` option for Python.  The Python module
will be placed in a subdirectory of your main build directory called “python”)

*Documentation for APIs in other languages:*  While there is extensive
Doxygen documentation available for the C++ and Python APIs, there is no
separate on-line documentation for the C and Fortran API. Instead, you should
use the C++ documentation, employing the mappings described here to figure out
the equivalent syntax in C or Fortran.

C API
*****

Before you start writing your own C program that calls OpenMM, be sure you can
build and run the two C examples that are supplied with OpenMM (see Chapter :numref:`openmm-tutorials`\ ).
These can be built from the supplied :code:`Makefile` on Linux and Mac, or
supplied :code:`NMakefile` and Visual Studio solution files on Windows.

The example programs are :code:`HelloArgonInC` and
:code:`HelloSodiumChlorideInC`\ . The argon example serves as a quick check that
your installation is set up properly and you know how to build a C program that
is linked with OpenMM. It will also tell you whether OpenMM is executing on the
GPU or is running (slowly) on the Reference platform. However, the argon example
is not a good template to follow for your own programs. The sodium chloride
example, though necessarily simplified, is structured roughly in the way we
recommended you set up your own programs to call OpenMM. Please be sure you have
both of these programs executing successfully on your machine before continuing.

Mechanics of using the C API
============================

The C API is generated automatically from the C++ API when OpenMM is built.
There are two resulting components: C bindings (functions to call), and C
declarations (in a header file). The C bindings are small :code:`extern`
(global) interface functions, one for every method of every OpenMM class, whose
signatures (name and arguments) are predictable from the class name and method
signatures. There are also “helper” types and functions provided for the few
cases in which the C++ behavior cannot be directly mapped into C. These
interface and helper functions are compiled in to the main OpenMM library so
there is nothing special you have to do to get access to them.

In the :code:`include` subdirectory of your OpenMM installation directory,
there is a machine-generated header file :code:`OpenMMCWrapper.h` that
should be #included in any C program that is to make calls to OpenMM functions.
That header contains declarations for all the OpenMM C interface functions and
related types. Note that if you follow our suggested structure, you will not
need to include this file in your :code:`main()` compilation unit but can
instead use it only in a local file that you write to provide a simple interface
to your existing code (see Chapter :numref:`openmm-tutorials`).

Mapping from the C++ API to the C API
=====================================

The automated generator of the C “wrappers” follows the translation strategy
shown in :autonumref:`Table,C API`\ . The idea is that if you see the construct on the left in
the C++ API documentation, you should interpret it as the corresponding
construct on the right in C. Please look at the supplied example programs to see
how this is done in practice.

==========================  =========================================  ===================================================
Construct                   C++ API declaration                        Equivalent in C API
==========================  =========================================  ===================================================
namespace                   OpenMM\::                                  OpenMM\_ (prefix)
class                       class OpenMM::ClassName                    typedef OpenMM_ClassName
constant                    OpenMM::RadiansPerDeg                      OpenMM_RadiansPerDeg (static constant)
class enum                  OpenMM::State::Positions                   OpenMM_State_Positions
constructor                 new OpenMM::ClassName()                    | OpenMM_ClassName* OpenMM_ClassName_create()
                                                                       | (additional constructors are _create_2(), etc.)
destructor                  | OpenMM::ClassName* thing;                | OpenMM_ClassName* thing;
                            | delete thing;                            | OpenMM_ClassName_destroy(thing);
class method                | OpenMM::ClassName* thing;                | OpenMM_ClassName* thing;
                            | thing->method(args);                     | OpenMM_ClassName_method(thing, args)
Boolean (type & constants)  | bool                                     | OpenMM_Boolean
                            | true, false                              | OpenMM_True(1), OpenMM_False(0)
string                      std::string                                char*
3-vector                    OpenMM::Vec3                               typedef OpenMM_Vec3
arrays                      | std::vector<std::string>                 | typedef OpenMM_StringArray
                            | std::vector<double>                      | typedef OpenMM_DoubleArray
                            | std::vector<Vec3>                        | typedef OpenMM_Vec3Array
                            | std::vector<std::pair<int,int>>          | typedef OpenMM_BondArray
                            | std::map<std::string,double>             | typedef OpenMM_ParameterArray
==========================  =========================================  ===================================================

:autonumber:`Table,C API`\ : Default mapping of objects from the C++ API to the C API
There are some exceptions to the generic translation rules shown in the table;
they are enumerated in the next section. And because there are no C++ API
equivalents to the array types, they are described in detail below.

Exceptions
==========

These two methods are handled somewhat differently in the C API than in the C++ API:

* **OpenMM::Context::getState()** The C version,
  :code:`OpenMM_Context_getState()`\ , returns a pointer to a heap allocated
  :code:`OpenMM_State` object. You must then explicitly destroy this
  :code:`State` object when you are done with it, by calling
  :code:`OpenMM_State_destroy()`\ .
* **OpenMM::Platform::loadPluginsFromDirectory()** The C version
  :code:`OpenMM_Platform_loadPluginsFromDirectory()` returns a heap-allocated
  :code:`OpenMM_StringArray` object containing a list of all the file names
  that were successfully loaded. You must then explicitly destroy this
  :code:`StringArray` object when you are done with it. Do not ignore the return
  value; if you do you’ll have a memory leak since the :code:`StringArray`
  will still be allocated.


(In the C++ API, the equivalent methods return references into existing memory
rather than new heap-allocated memory, so the returned objects do not need to be
destroyed.)

OpenMM_Vec3 helper type
=======================

Unlike the other OpenMM objects which are opaque and manipulated via pointers,
the C API provides an explicit definition for the C :code:`OpenMM_Vec3` type
that is compatible with the :code:`OpenMM::Vec3` type. The definition of
:code:`OpenMM_Vec3` is:

.. code-block:: c

    typedef struct {double x, y, z;} OpenMM_Vec3;

You can work directly with the individual fields of this type from your C
program if you want. For convenience, a scale() function is provided that
creates a new OpenMM_Vec3 from an old one and a scale factor:

.. code-block:: c

    OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale);

Array helper types
==================

C++ has built-in container types :code:`std::vector` and :code:`std::map`
which OpenMM uses to manipulate arrays of objects. These don’t have direct
equivalents in C, so we supply special array types for each kind of object for
which OpenMM creates containers. These are: string, double, Vec3, bond, and
parameter map. See :autonumref:`Table,C arrays` for the names of the C types for each of these
object arrays. Each of the array types provides these functions (prefixed by
:code:`OpenMM_` and the actual *Thing* name), with the syntax shown
conceptually since it differs slightly for each kind of object.

.. tabularcolumns:: |l|L|

=======================================================  =========================================================================================================================================================================================================
Function                                                 Operation
=======================================================  =========================================================================================================================================================================================================
*Thing*\ Array\* create(int size)                        Create a heap-allocated array of *Things*\ , with space pre-allocated to hold :code:`size` of them. You can start at :code:`size==0` if you want since these arrays are dynamically resizeable.
void destroy(\ *Thing*\ Array\*)                         Free the heap space that is currently in use for the passed-in array of *Things*\ .
int getSize(\ *Thing*\ Array\*)                          Return the current number of *Things* in this array. This means you can :code:`get()` and :code:`set()` elements up to :code:`getSize()-1`\ .
void resize(\ *Thing*\ Array\*, int size)                Change the size of this array to the indicated value which may be smaller or larger than the current size. Existing elements remain in their same locations as long as they still fit.
void append(\ *Thing*\ Array\*, *Thing*\ )               Add a *Thing* to the end of the array, increasing the array size by one. The precise syntax depends on the actual type of *Thing*\ ; see below.
void set(\ *Thing*\ Array\*, int index, *Thing*\ )       Store a copy of *Thing* in the indicated element of the array (indexed from 0). The array must be of length at least :code:`index+1`\ ; you can’t grow the array with this function.
*Thing* get(\ *Thing*\ Array\*, int index)               Retrieve a particular element from the array (indexed from 0). (For some Things the value is returned in arguments rather than as the function return.)
=======================================================  =========================================================================================================================================================================================================

:autonumber:`Table,C arrays`\ : Generic description of array helper types

Here are the exact declarations with deviations from the generic description
noted, for each of the array types.

OpenMM_DoubleArray
------------------

.. code-block:: c

    OpenMM_DoubleArray*
                OpenMM_DoubleArray_create(int size);
    void        OpenMM_DoubleArray_destroy(OpenMM_DoubleArray*);
    int         OpenMM_DoubleArray_getSize(const OpenMM_DoubleArray*);
    void        OpenMM_DoubleArray_resize(OpenMM_DoubleArray*, int size);
    void        OpenMM_DoubleArray_append(OpenMM_DoubleArray*, double value);
    void        OpenMM_DoubleArray_set(OpenMM_DoubleArray*, int index, double value);
    double      OpenMM_DoubleArray_get(const OpenMM_DoubleArray*, int index);

OpenMM_StringArray
------------------

.. code-block:: c

    OpenMM_StringArray*
                OpenMM_StringArray_create(int size);
    void        OpenMM_StringArray_destroy(OpenMM_StringArray*);
    int         OpenMM_StringArray_getSize(const OpenMM_StringArray*);
    void        OpenMM_StringArray_resize(OpenMM_StringArray*, int size);
    void        OpenMM_StringArray_append(OpenMM_StringArray*, const char* string);
    void        OpenMM_StringArray_set(OpenMM_StringArray*, int index, const char* string);
    const char* OpenMM_StringArray_get(const OpenMM_StringArray*, int index);

OpenMM_Vec3Array
----------------

.. code-block:: c

    OpenMM_Vec3Array*
                OpenMM_Vec3Array_create(int size);
    void        OpenMM_Vec3Array_destroy(OpenMM_Vec3Array*);
    int         OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array*);
    void        OpenMM_Vec3Array_resize(OpenMM_Vec3Array*, int size);
    void        OpenMM_Vec3Array_append(OpenMM_Vec3Array*, const OpenMM_Vec3 vec);
    void        OpenMM_Vec3Array_set(OpenMM_Vec3Array*, int index, const OpenMM_Vec3 vec);
    const OpenMM_Vec3*
                OpenMM_Vec3Array_get(const OpenMM_Vec3Array*, int index);

OpenMM_BondArray
----------------

Note that bonds are specified by pairs of integers (the atom indices). The
:code:`get()` method returns those in a pair of final arguments rather than as
its functional return.

.. code-block:: c

    OpenMM_BondArray*
                OpenMM_BondArray_create(int size);
    void        OpenMM_BondArray_destroy(OpenMM_BondArray*);
    int         OpenMM_BondArray_getSize(const OpenMM_BondArray*);
    void        OpenMM_BondArray_resize(OpenMM_BondArray*, int size);
    void        OpenMM_BondArray_append(OpenMM_BondArray*, int particle1, int particle2);
    void        OpenMM_BondArray_set(OpenMM_BondArray*, int index, int particle1, int particle2);
    void        OpenMM_BondArray_get(const OpenMM_BondArray*, int index,
                                     int* particle1, int* particle2);

OpenMM_ParameterArray
---------------------

OpenMM returns references to internal :code:`ParameterArrays` but does not
support user-created :code:`ParameterArrays`\ , so only the :code:`get()`
and :code:`getSize()` functions are available. Also, note that since this is
actually a map rather than an array, the “index” is the *name* of the
parameter rather than its ordinal.

.. code-block:: c

    int         OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray*);
    double      OpenMM_ParameterArray_get(const OpenMM_ParameterArray*, const char* name);


Fortran 95 API
*****************

Before you start writing your own Fortran program that calls OpenMM, be sure you
can build and run the two Fortran examples that are supplied with OpenMM (see
Chapter :numref:`openmm-tutorials`). These can be built from the supplied :code:`Makefile` on Linux
and Mac, or supplied :code:`NMakefile` and Visual Studio solution files on
Windows.

The example programs are :code:`HelloArgonInFortran` and
:code:`HelloSodiumChlorideInFortran`\ . The argon example serves as a quick
check that your installation is set up properly and you know how to build a
Fortran program that is linked with OpenMM. It will also tell you whether OpenMM
is executing on the GPU or is running (slowly) on the Reference platform.
However, the argon example is not a good template to follow for your own
programs. The sodium chloride example, though necessarily simplified, is
structured roughly in the way we recommended you set up your own programs to
call OpenMM. Please be sure you have both of these programs executing
successfully on your machine before continuing.

Mechanics of using the Fortran API
==================================

The Fortran API is generated automatically from the C++ API when OpenMM is
built. There are two resulting components: Fortran bindings (subroutines to
call), and Fortran declarations of types and subroutines (in the form of a
Fortran 95 module file). The Fortran bindings are small interface subroutines,
one for every method of every OpenMM class, whose signatures (name and
arguments) are predictable from the class name and method signatures. There are
also “helper” types and subroutines provided for the few cases in which the C++
behavior cannot be directly mapped into Fortran. These interface and helper
subroutines are compiled in to the main OpenMM library so there is nothing
special you have to do to get access to them.

Because Fortran is case-insensitive, calls to Fortran subroutines (however
capitalized) are mapped by the compiler into all-lowercase or all-uppercase
names, and different compilers use different conventions. The automatically-generated
OpenMM Fortran “wrapper” subroutines, which are generated in C and
thus case-sensitive, are provided in two forms for compatibility with the
majority of Fortran compilers, including Intel Fortran and gfortran. The two
forms are: (1) all-lowercase with a trailing underscore, and (2) all-uppercase
without a trailing underscore. So regardless of the Fortran compiler you are
using, it should find a suitable subroutine to call in the main OpenMM library.

In the :code:`include` subdirectory of your OpenMM installation directory,
there is a machine-generated module file :code:`OpenMMFortranModule.f90`
that must be compiled along with any Fortran program that is to make calls to
OpenMM functions. (You can look at the :code:`Makefile` or Visual Studio
solution file provided with the OpenMM examples to see how to build a program
that uses this module file.) This module file contains definitions for two
modules: :code:`MODULE OpenMM_Types` and :code:`MODULE OpenMM`\ ; however,
only the :code:`OpenMM` module will appear in user programs (it references
the other module internally). The modules contain declarations for all the
OpenMM Fortran interface subroutines, related types, and parameters (constants).
Note that if you follow our suggested structure, you will not need to
:code:`use` the :code:`OpenMM` module in your :code:`main()`
compilation unit but can instead use it only in a local file that you write to
provide a simple interface to your existing code (see Chapter :numref:`openmm-tutorials`).

Mapping from the C++ API to the Fortran API
===========================================

The automated generator of the Fortran “wrappers” follows the translation
strategy shown in :autonumref:`Table,Fortran API`\ . The idea is that if you see the construct on the
left in the C++ API documentation, you should interpret it as the corresponding
construct on the right in Fortran. Please look at the supplied example programs
to see how this is done in practice. Note that all subroutines and modules are
declared with “\ :code:`implicit none`\ ”, meaning that the type of every symbol
is declared explicitly and should not be inferred from the first letter of the
symbol name.

==========================  ===================================  ========================================================
Construct                   C++ API declaration                  Equivalent in Fortran API
==========================  ===================================  ========================================================
namespace                   OpenMM\::                            OpenMM\_ (prefix)
class                       class OpenMM::ClassName              type (OpenMM_ClassName)
constant                    OpenMM::RadiansPerDeg                parameter (OpenMM_RadiansPerDeg)
class enum                  OpenMM::State::Positions             parameter (OpenMM_State_Positions)
constructor                 new OpenMM::ClassName()              | type (OpenMM_ClassName) thing
                                                                 | call OpenMM_ClassName_create(thing)
                                                                 | (additional constructors are \_create_2(), etc.)
destructor                  | OpenMM::ClassName* thing;          | type (OpenMM_ClassName) thing
                            | delete thing;                      | call OpenMM_ClassName_destroy(thing)
class method                | OpenMM::ClassName* thing;          | type (OpenMM_ClassName) thing
                            | thing->method(args*)               | call OpenMM_ClassName_method(thing, args)
Boolean (type & constants)  | bool                               | integer*4
                            | true                               | parameter (OpenMM_True=1)
                            | false                              | parameter (OpenMM_False=0)
string                      std::string                          character(*)
3-vector                    OpenMM::Vec3                         real*8 vec(3)
arrays                      std::vector<std::string>             | type (OpenMM_StringArray)
                            std::vector<double>                  | type (OpenMM_DoubleArray)
                            std::vector<Vec3>                    | type (OpenMM_Vec3Array)
                            std::vector<std::pair<int,int>>      | type (OpenMM_BondArray)
                            std::map<std::string, double>        | type (OpenMM_ParameterArray)
==========================  ===================================  ========================================================

:autonumber:`Table,Fortran API`\ : Default mapping of objects from the C++ API to the Fortran API

Because there are no C++ API equivalents to the array types, they are described
in detail below.

OpenMM_Vec3 helper type
=======================

Unlike the other OpenMM objects which are opaque and manipulated via pointers,
the Fortran API uses an ordinary :code:`real*8(3)` array in
place of the :code:`OpenMM::Vec3` type.
You can work directly with the individual elements of this type from your
Fortran program if you want. For convenience, a :code:`scale()` function is
provided that creates a new Vec3 from an old one and a scale factor:

.. code-block:: fortran

    subroutine OpenMM_Vec3_scale(vec, scale, result)
    real*8 vec(3), scale, result(3)

No explicit :code:`type`\ :code:`(OpenMM_Vec3)` is provided in the Fortran
API since it is not needed.

Array helper types
==================

C++ has built-in container types :code:`std::vector` and :code:`std::map`
which OpenMM uses to manipulate arrays of objects. These don’t have direct
equivalents in Fortran, so we supply special array types for each kind of object
for which OpenMM creates containers. These are: string, double, Vec3, bond, and
parameter map. See :autonumref:`Table,Fortran arrays` for the names of the Fortran types for each of
these object arrays. Each of the array types provides these functions (prefixed
by :code:`OpenMM_` and the actual *Thing* name), with the syntax shown
conceptually since it differs slightly for each kind of object.

+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| Function                                  | Operation                                                                                              |
+===========================================+========================================================================================================+
| | subroutine create(array,size)           | Create a heap-allocated array of *Things*\ , with space pre-allocated to hold :code:`size` of them.    |
| | type (OpenMM\_\ *Thing*\ Array) array   | You can start at :code:`size`\ ==0 if you want since these arrays are dynamically resizeable.          |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine destroy(array)               | Free the heap space that is currently in use for the passed-in array of *Things*\ .                    |
| | type (OpenMM\_\ *Thing*\ Array) array   |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | function getSize(array)                 | Return the current number of *Things* in this array. This means you can :code:`get()` and              |
| | type (OpenMM\_\ *Thing*\ Array) array   | :code:`set()` elements up to :code:`getSize()`\ .                                                      |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine resize(array,size)           | Change the size of this array to the indicated value which may be smaller or larger than the           |
| | type (OpenMM\_\ *Thing*\ Array) array   | current size. Existing elements remain in their same locations as long as they still fit.              |
| | integer*4 size                          |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine append(array,elt)            | Add a *Thing* to the end of the array, increasing the array size by one. The precise syntax depends    |
| | type (OpenMM\_\ *Thing*\ Array) array   | on the actual type of *Thing*\ ; see below.                                                            |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine set(array,index,elt)         | Store a copy of :code:`elt` in the indicated element of the array (indexed from 1). The array must     |
| | type (OpenMM\_\ *Thing*\ Array) array   | be of length at least :code:`index`\ ; you can’t grow the array with this function.                    |
| | integer*4 size                          |                                                                                                        |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+
| | subroutine get(array,index,elt)         | Retrieve a particular element from the array (indexed from 1).  Some *Things* require more than one    |
| | type (OpenMM\_\ *Thing*\ Array) array   | argument to return.                                                                                    |
| | integer*4 size                          |                                                                                                        |
| | *Thing* elt                             |                                                                                                        |
+-------------------------------------------+--------------------------------------------------------------------------------------------------------+

:autonumber:`Table,Fortran arrays`\ : Generic description of array helper types

Here are the exact declarations with deviations from the generic description
noted, for each of the array types.

OpenMM_DoubleArray
------------------

.. code-block:: fortran

    subroutine OpenMM_DoubleArray_create(array, size)
        integer*4 size
        type (OpenMM_DoubleArray) array
    subroutine OpenMM_DoubleArray_destroy(array)
        type (OpenMM_DoubleArray) array
    function OpenMM_DoubleArray_getSize(array)
        type (OpenMM_DoubleArray) array
        integer*4 OpenMM_DoubleArray_getSize
    subroutine OpenMM_DoubleArray_resize(array, size)
        type (OpenMM_DoubleArray) array
        integer*4 size
    subroutine OpenMM_DoubleArray_append(array, value)
        type (OpenMM_DoubleArray) array
        real*8 value
    subroutine OpenMM_DoubleArray_set(array, index, value)
        type (OpenMM_DoubleArray) array
        integer*4 index
        real*8 value
    subroutine OpenMM_DoubleArray_get(array, index, value)
        type (OpenMM_DoubleArray) array
        integer*4 index
        real*8 value

OpenMM_StringArray
------------------

.. code-block:: fortran

    subroutine OpenMM_StringArray_create(array, size)
        integer*4 size
        type (OpenMM_StringArray) array
    subroutine OpenMM_StringArray_destroy(array)
        type (OpenMM_StringArray) array
    function OpenMM_StringArray_getSize(array)
        type (OpenMM_StringArray) array
        integer*4 OpenMM_StringArray_getSize
    subroutine OpenMM_StringArray_resize(array, size)
        type (OpenMM_StringArray) array
        integer*4 size
    subroutine OpenMM_StringArray_append(array, str)
        type (OpenMM_StringArray) array
        character(*) str
    subroutine OpenMM_StringArray_set(array, index, str)
        type (OpenMM_StringArray) array
        integer*4 index
        character(*) str
    subroutine OpenMM_StringArray_get(array, index, str)
        type (OpenMM_StringArray) array
        integer*4 index
        character(*)str

OpenMM_Vec3Array
----------------

.. code-block:: fortran

    subroutine OpenMM_Vec3Array_create(array, size)
        integer*4 size
        type (OpenMM_Vec3Array) array
    subroutine OpenMM_Vec3Array_destroy(array)
        type (OpenMM_Vec3Array) array
    function OpenMM_Vec3Array_getSize(array)
        type (OpenMM_Vec3Array) array
        integer*4 OpenMM_Vec3Array_getSize
    subroutine OpenMM_Vec3Array_resize(array, size)
        type (OpenMM_Vec3Array) array
        integer*4 size
    subroutine OpenMM_Vec3Array_append(array, vec)
        type (OpenMM_Vec3Array) array
        real*8 vec(3)
    subroutine OpenMM_Vec3Array_set(array, index, vec)
        type (OpenMM_Vec3Array) array
        integer*4 index
        real*8 vec(3)
    subroutine OpenMM_Vec3Array_get(array, index, vec)
        type (OpenMM_Vec3Array) array
        integer*4 index
        real*8 vec (3)

OpenMM_BondArray
----------------

Note that bonds are specified by pairs of integers (the atom indices). The
:code:`get()` method returns those in a pair of final arguments rather than as
its functional return.

.. code-block:: fortran

    subroutine OpenMM_BondArray_create(array, size)
        integer*4 size
        type (OpenMM_BondArray) array
    subroutine OpenMM_BondArray_destroy(array)
        type (OpenMM_BondArray) array
    function OpenMM_BondArray_getSize(array)
        type (OpenMM_BondArray) array
        integer*4 OpenMM_BondArray_getSize
    subroutine OpenMM_BondArray_resize(array, size)
        type (OpenMM_BondArray) array
        integer*4 size
    subroutine OpenMM_BondArray_append(array, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 particle1, particle2
    subroutine OpenMM_BondArray_set(array, index, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 index, particle1, particle2
    subroutine OpenMM_BondArray_get(array, index, particle1, particle2)
        type (OpenMM_BondArray) array
        integer*4 index, particle1, particle2

OpenMM_ParameterArray
---------------------

OpenMM returns references to internal :code:`ParameterArrays` but does not
support user-created :code:`ParameterArrays`\ , so only the :code:`get()`
and :code:`getSize()` functions are available. Also, note that since this is
actually a map rather than an array, the “index” is the *name* of the
parameter rather than its ordinal.

.. code-block:: fortran

    function OpenMM_ParameterArray_getSize(array)
        type (OpenMM_ParameterArray) array
        integer*4 OpenMM_ParameterArray_getSize
    subroutine OpenMM_ParameterArray_get(array, name, param)
        type (OpenMM_ParameterArray) array
        character(*) name
        character(*) param


Python API
**********


Mapping from the C++ API to the Python API
==========================================

The Python API follows the C++ API as closely as possible. There are three
notable differences:

#. The :code:`getState()` method in the :code:`Context` class takes
   Pythonic-type arguments to indicate which state variables should be made
   available.  For example:
   ::

    myContext.getState(energy=True, force=False, …)

#. Wherever the C++ API uses references to return multiple values from a method,
   the Python API returns a tuple.  For example, in C++ you would query a
   HarmonicBondForce for a bond’s parameters as follows:
   ::

    int particle1, particle2;
    double length, k;
    f.getBondParameters(i, particle1, particle2, length, k);

   In Python, the equivalent code is:
   ::

    [particle1, particle2, length, k] = f.getBondParameters(i)

#. Unlike C++, the Python API accepts and returns quantities with units attached
   to most values (see Section :numref:`units-and-dimensional-analysis` below for
   details).  In short, this means that while values in C++ have *implicit*
   units, the Python API returns objects that have values and *explicit* units.


Mechanics of using the Python API
=================================

When using the Python API, be sure to include the GPU support
libraries in your library path, just as you would for a C++ application.  This
is set with the :code:`LD_LIBRARY_PATH` environment variable on Linux,
:code:`DYLD_LIBRARY_PATH` on Mac, or :code:`PATH` on Windows.  See
Chapter :numref:`installing-openmm` for details.

The Python API is contained in the openmm package, while the units code is
contained in the openmm.units package.  (The application layer, described in the
Application Guide, is contained in the openmm.app package.)  A program
using it will therefore typically begin
::

    import openmm as mm
    import openmm.unit as unit

Creating and using OpenMM objects is then done exactly as in C++:
::

    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(1.2*unit.nanometer)
    system.addForce(nb)

Note that when setting the cutoff distance, we explicitly specify that it is in
nanometers.  We could just as easily specify it in different units:
::

    nb.setCutoffDistance(12*unit.angstrom)

The use of units in OpenMM is discussed in the next section.


.. _units-and-dimensional-analysis:

Units and dimensional analysis
==============================


Why does the Python API include units?
--------------------------------------

The C++ API for OpenMM uses an *implicit* set of units for physical
quantities such as lengths, masses, energies, etc.  These units are based on
daltons, nanometers, and picoseconds for the mass, length, and time dimensions,
respectively.  When using the C++ API, it is very important to ensure that
quantities being manipulated are always expressed in terms of these units.  For
example, if you read in a distance in Angstroms, you must multiply that distance
by a conversion factor to turn it into nanometers before using it in the C++
API.  Such conversions can be a source of tedium and errors.  This is true in
many areas of scientific programming.  Units confusion was blamed for the loss
of the Mars Climate Orbiter spacecraft in 1999, at a cost of more than $100
million.  Units were introduced in the Python API to minimize the chance of such
errors.

The Python API addresses the potential problem of conversion errors by using
quantities with explicit units.  If a particular distance is expressed in
Angstroms, the Python API will know that it is in Angstroms.  When the time
comes to call the C++ API, it will understand that the quantity must be
converted to nanometers.  You, the programmer, must declare upfront that the
quantity is in Angstrom units, and the API will take care of the details from
then on.  Using explicit units is a bit like brushing your teeth: it requires
some effort upfront, but it probably saves you trouble in the long run.

Quantities, units, and dimensions
---------------------------------

The explicit unit system is based on three concepts: Dimensions, Units, and
Quantities.

Dimensions are measurable physical concepts such as mass, length, time, and
energy.  Energy is actually a composite dimension based on mass, length, and
time.

A Unit defines a linear scale used to measure amounts of a particular physical
Dimension.  Examples of units include meters, seconds, joules, inches, and
grams.

A Quantity is a specific amount of a physical Dimension.  An example of a
quantity is “0.63 kilograms”.  A Quantity is expressed as a combination of a
value (e.g., 0.63), and a Unit (e.g., kilogram).  The same Quantity can be
expressed in different Units.

The set of BaseDimensions defined in the openmm.unit module includes:

* mass
* length
* time
* temperature
* amount
* charge
* luminous intensity


These are not precisely the same list of base dimensions used in the SI unit
system.  SI defines “current” (charge per time) as a base unit, while openmm.unit
uses “charge”.  And openmm.unit treats angle as a dimension, even though angle
quantities are often considered dimensionless.  In this case, we choose to err
on the side of explicitness, particularly because interconversion of degrees and
radians is a frequent source of unit headaches.

Units examples
--------------

Many common units are defined in the openmm.unit module.
::

    from openmm.unit import nanometer, angstrom, dalton

Sometimes you don’t want to type the full unit name every time, so you can
assign it a shorter name using the :code:`as` functionality:
::

    from openmm.unit import nanometer as nm

New quantities can be created from a value and a unit.  You can use either the
multiply operator (‘*’) or the explicit Quantity constructor:
::

    from simk.unit import nanometer, Quantity
    # construct a Quantity using the multiply operator
    bond_length = 1.53 * nanometer
    # equivalently using the explicit Quantity constructor
    bond_length = Quantity(1.53, nanometer)
    # or more verbosely
    bond_length = Quantity(value=1.53, unit=nanometer)

Arithmetic with units
---------------------

Addition and subtraction of quantities is only permitted between quantities that
share the same dimension.  It makes no sense to add a mass to a distance.  If
you attempt to add or subtract two quantities with different dimensions, an
exception will be raised.  This is a good thing; it helps you avoid errors.
::

    x = 5.0*dalton + 4.3*nanometer; # error

Addition or subtraction of quantities with the same dimension, but different
units, is fine, and results in a new quantity created using the correct
conversion factor between the units used.
::

    x = 1.3*nanometer + 5.6*angstrom; # OK, result in nanometers

Quantities can be added and subtracted.  Naked Units cannot.

Multiplying or dividing two quantities creates a new quantity with a composite
dimension.  For example, dividing a distance by a time results in a velocity.
::

    from openmm.unit import kilogram, meter, second
    a = 9.8 * meter / second**2; # acceleration
    m = 0.36 * kilogram; # mass
    F = m * a; # force in kg*m/s**2::


Multiplication or division of two Units results in a composite Unit.
::

    mps = meter / second

Unlike amount (moles), angle (radians) is arguably dimensionless.  But openmm.unit
treats angle as another dimension.   Use the trigonometric functions from the
openmm.unit module (not those from the Python math module!) when dealing with
Units and Quantities.
::

    from openmm.unit import sin, cos, acos
    x = sin(90.0*degrees)
    angle = acos(0.68); # returns an angle quantity (in radians)

The method :code:`pow()` is a built-in Python method that works with
Quantities and Units.
::

    area = pow(3.0*meter, 2)
    # or, equivalently
    area = (3.0*meter)**2
    # or
    area = 9.0*(meter**2)

The method :code:`sqrt()` is not as built-in as :code:`pow()`\ .  Do not
use the Python :code:`math.sqrt()` method with Units and Quantities.  Use
the :code:`openmm.unit.sqrt()` method instead:
::

    from openmm.unit import sqrt
    side_length = sqrt(4.0*meter**2)


Atomic scale mass and energy units are “per amount”
---------------------------------------------------

Mass and energy units at the atomic scale are specified “per amount” in the
openmm.unit module.  Amount (mole) is one of the seven fundamental dimensions in
the SI unit system.   The atomic scale mass unit, dalton, is defined as grams
per mole.  The dimension of dalton is therefore mass/amount, instead of simply
mass.  Similarly, the atomic scale energy unit, kilojoule_per_mole (and
kilocalorie_per_mole) has “per amount” in its dimension.  Be careful to always
use “per amount” mass and energy types at the atomic scale, and your dimensional
analysis should work out properly.

The energy unit kilocalories_per_mole does not have the same Dimension as the
macroscopic energy unit kilocalories.  Molecular scientists sometimes use the
word "kilocalories" when they mean "kilocalories per mole".  Use "kilocalories
per mole" or"kilojoules per mole" for molecular energies.  Use "kilocalories"
for the metabolic energy content of your lunch.  The energy unit
kilojoule_per_mole happens to go naturally with the units nanometer,
picoseconds, and dalton.  This is because 1 kilojoule/mole happens to be equal
to 1 gram-nanometer\ :sup:`2`\ /mole-picosecond\ :sup:`2`\ , and is therefore
consistent with the molecular dynamics unit system used in the C++ OpenMM API.

These "per mole" units are what you should be using for molecular calculations,
as long as you are using SI / cgs / calorie sorts of units.

SI prefixes
-----------

Many units with SI prefixes such as “milligram” (milli) and “kilometer” (kilo)
are provided in the openmm.unit module.  Others can be created by multiplying a
prefix symbol by a non-prefixed unit:
::

    from openmm.unit import mega, kelvin
    megakelvin = mega * kelvin
    t = 8.3 * megakelvin

Only grams and meters get all of the SI prefixes (from yotto-(10\ :sup:`-24`\ )
to yotta-(10\ :sup:`24`\ )) automatically.


Converting to different units
-----------------------------

Use the :code:`Quantity.in_units_of()` method to create a new Quantity with
different units.
::

    from openmm.unit import nanosecond, fortnight
    x = (175000*nanosecond).in_units_of(fortnight)

When you want a plain number out of a Quantity, use the :code:`value_in_unit()` method:
::

    from openmm.unit import femtosecond, picosecond
    t = 5.0*femtosecond
    t_just_a_number = t.value_in_unit(picoseconds)

Using :code:`value_in_unit()` puts the responsibility for unit analysis back
into your hands, and it should be avoided.  It is sometimes necessary, however,
when you are called upon to use a non-units-aware Python API.


Lists, tuples, vectors, numpy arrays, and Units
-----------------------------------------------

Units can be attached to containers of numbers to create a vector quantity.  The
openmm.unit module overloads the :code:`__setitem__` and
:code:`__getitem__` methods for these containers to ensure that Quantities go
in and out.
::

    >>> a = Vec3(1,2,3) * nanometers
    >>> print(a)
    (1, 2, 3) nm
    >>> print(a.in_units_of(angstroms))
    (10.0, 20.0, 30.0) A

    >>> s2 = [[1,2,3],[4,5,6]] * centimeter
    >>> print(s2)
    [[1, 2, 3], [4, 5, 6]] cm
    >>> print(s2/millimeter)
    [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]

    >>> import numpy
    >>> a = numpy.array([1,2,3]) * centimeter
    >>> print(a)
    [1 2 3] cm
    >>> print(a/millimeter)
    [ 10.  20.  30.]

Converting a whole list to different units at once is much faster than
converting each element individually.  For example, consider the following code
that prints out the position of every particle in a State, as measured in
Angstroms:
::

    for v in state.getPositions():
        print(v.value_in_unit(angstrom))

This can be rewritten as follows:
::

    for v in state.getPositions().value_in_unit(angstrom):
        print(v)

The two versions produce identical results, but the second one will run faster,
and therefore is preferred.
