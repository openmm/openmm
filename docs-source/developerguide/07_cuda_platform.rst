.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _the-cuda-platform:

The CUDA and HIP Platforms
##########################

The CUDA and HIP platforms are very similar to the OpenCL platform, and most of the
previous chapter applies equally well to them, just changing “OpenCL” to “Cuda” or
"Hip" in class names.  There are a few differences worth noting.

Caching Kernels
***************

Like the OpenCL platform, the CUDA and HIP platforms compile their kernels at runtime.
To improve performance, they try to cache the compiled kernels on disk for
later use.  This allows subsequent Contexts to skip compiling some kernels.  To
make this work, they need a directory on disk where they can write out temporary
files.  It is specified by the “TempDirectory” property when you create a
new Context.  They usually can figure out a suitable value on their own, but
sometimes they need help.  See the “Platform-Specific Properties” chapter of the
User's Manual for details.

Accumulating Forces
*******************

The OpenCL platform, as described in Section :numref:`computing-forces`\ , uses two types of buffers for
accumulating forces: a set of floating point buffers, and a single fixed point
buffer.  In contrast, the CUDA and HIP platforms use *only* the fixed point buffer
(represented by the CUDA type :code:`long` :code:`long`\ ).
