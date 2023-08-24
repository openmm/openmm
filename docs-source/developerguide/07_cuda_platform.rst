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

Caching Kernels
***************

Like the OpenCL platform, the CUDA platform compiles all its kernels at runtime.
To improve performance, it tries to cache the compiled kernels on disk for
later use.  This allows subsequent Contexts to skip compiling some kernels.  To
make this work, it needs a directory on disk where it can write out temporary
files.  It is specified by the “CudaTempDirectory” property when you create a
new Context.  It usually can figure out a suitable value on its own, but
sometimes it needs help.  See the “Platform-Specific Properties” chapter of the
User's Manual for details.

Accumulating Forces
*******************

The OpenCL platform, as described in Section :numref:`computing-forces`\ , uses two types of buffers for
accumulating forces: a set of floating point buffers, and a single fixed point
buffer.  In contrast, the CUDA platform uses *only* the fixed point buffer
(represented by the CUDA type :code:`long` :code:`long`\ ).
