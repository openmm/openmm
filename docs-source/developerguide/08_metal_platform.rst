.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-metal-platform:

The Metal Platform
##################

The Metal platform is very similar to the CUDA and HIP platforms.  Due to particular
features of Apple GPUs and the Metal framework, there are a few differences to be
aware of.

Single GPU Only
***************

The Metal platform only runs on Macintosh computers with ARM processors.  None of those computers
supports multiple GPUs.  The Metal platform therefore does not have an option to select a GPU, and
does not include the code for parallelizing simulations across multiple GPUs found in the other
platforms.

Passing Arguments to Kernels
****************************

Metal has restrictions on the arguments passed to kernels that are not found in any other framework.
First, it does not allow arguments to be primitive types like int or float.  Any primitive type
argument must be transformed into a reference to a value in constant memory.  Second, it only allows
a maximum of 31 arguments to any kernel.  When more arguments are required, they must be packed into
a larger object.

OpenMM's Metal platform tries to handle these complexities automatically.  When you give it the
source code for a program to compile, it looks for kernel signatures and automatically rewrites
them to work around the problems.  Usually this works perfectly and you do not need to think about
it.  It is good to be aware of it, though, in case you do run into a problem.

Mixed Precision
***************

Apple GPUs do not support double precision math.  OpenMM's Metal platform therefore provides an
implementation of the DF64 number format, which uses a pair of single precision numbers to
represent a single value at close to double precision.  This allows it to support mixed precision
mode despite the lack of native double precision.

DF64 largely functions as a drop in replacement for double precision, but be aware that it is
not exactly the same.  Its exact precision depends on the value being stored, but is always
at least 48 bits.  That is compared to 53 bits for double precision.  Also, its range is the
same as single precision (a maximum value of approximately :math:`3.4 \cdot 10^{38}`), which is
much smaller than the maximum double precision value.

Conversions between double and DF64 happen when setting kernel arguments, and when uploading or
downloading arrays.  There are a few things you need to do to let this happen correctly.  First,
when calling :code:`initialize()` on a ComputeArray object that holds double precision values,
be sure to call the version that takes a type as a template argument.  This allows the array to
recognize and record what type it contains.  Second, when calling :code:`upload()` or
:code:`download()` on the array, pass :code:`true` for the :code:`convert` argument.  This
instructs it to perform conversions between different floating point formats as needed.