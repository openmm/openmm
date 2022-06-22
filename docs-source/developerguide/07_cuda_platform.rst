.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-cuda-platform:

The CUDA Platform
#################

The CUDA platform is very similar to the OpenCL platform, and most of the
previous chapter applies equally well to it, just changing “OpenCL” to “Cuda” in
class names.  There are a few differences worth noting.

Compiling Kernels
*****************

Like the OpenCL platform, the CUDA platform compiles all its kernels at runtime.
Unlike OpenCL, CUDA does not have built in support for runtime compilation.
OpenMM therefore needs to implement this itself by writing the source code out
to disk, invoking the nvcc compiler as a separate process, and then loading the
compiled kernel in from disk.

For the most part, you can ignore all of this.  Just call
:code:`createModule()` on the CudaContext, passing it the CUDA source code.
It takes care of the details of compilation and loading, returning a CUmodule
object when it is done.  You can then call :code:`getKernel()` to look up
individual kernels in the module (represented as CUfunction objects) and
:code:`executeKernel()` to execute them.

The CUDA platform does need two things to make this work: a directory on disk
where it can write out temporary files, and the path to the nvcc compiler.
These are specified by the “CudaTempDirectory” and “CudaCompiler” properties
when you create a new Context.  It often can figure out suitable values for them
on its own, but sometimes it needs help.  See the “Platform-Specific Properties”
chapter of the User's Manual for details.

Accumulating Forces
*******************

The OpenCL platform, as described in Section :numref:`computing-forces`\ , uses two types of buffers for
accumulating forces: a set of floating point buffers, and a single fixed point
buffer.  In contrast, the CUDA platform uses *only* the fixed point buffer
(represented by the CUDA type :code:`long` :code:`long`\ ).  This means
the CUDA platform only works on devices that support 64 bit atomic operations
(compute capability 1.2 or higher).
