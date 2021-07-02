.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _common-compute:

Common Compute
##############

Common Compute is not a platform, but it shares many elements of one.  It exists
to reduce code duplication between the OpenCL and CUDA platforms.  It allows a
single implementation to be written for most kernels that can be used by both
platforms.

OpenCL and CUDA are very similar to each other.  Their computational models are
nearly identical.  For example, each is based around launching kernels that are
executed in parallel by many threads.  Each of them groups threads into blocks,
with more communication and synchronization permitted between the threads
in a block than between ones in different blocks.  They have very similar memory
hierarchies: high latency global memory, low latency local/shared memory that
can be used for communication between the threads of a block, and local variables
that are visible only to a single thread.

Even their languages for writing kernels are very similar.  Here is an OpenCL
kernel that adds two arrays together, storing the result in a third array.
::

    __kernel void addArrays(__global const float* restrict a,
                            __global const float* restrict b,
                            __global float* restrict c
                            int length) {
        for (int i = get_global_id(0); i < length; i += get_global_size(0))
            c[i] = a[i]+b[i];
    }

Here is the corresponding CUDA kernel.
::

    __extern "C" __global__ void addArrays(const float* __restrict__ a,
                                           const float* __restrict__ b,
                                           _float* __restrict__ c
                                           int length) {
        for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < length; i += blockDim.x*gridDim.x)
            c[i] = a[i]+b[i];
    }

The difference between them is largely just a mechanical find-and-replace.
After many years of writing and maintaining nearly identical kernels by hand,
it finally occurred to us that the translation could be done automatically by
the compiler.  Simply by defining a few preprocessor macros, the following
kernel can be compiled equally well either as OpenCL or as CUDA.
::

    KERNEL void addArrays(GLOBAL const float* RESTRICT a,
                          GLOBAL const float* RESTRICT b,
                          GLOBAL float* RESTRICT c
                          int length) {
        for (int i = GLOBAL_ID; i < length; i += GLOBAL_SIZE)
            c[i] = a[i]+b[i];
    }

Writing Device Code
*******************

When compiling kernels with the Common Compute API, the following macros are
defined.

+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|Macro                          |OpenCL Definition                                           |CUDA Definition                             |
+===============================+============================================================+============================================+
|:code:`KERNEL`                 |:code:`__kernel`                                            |:code:`extern "C" __global__`               |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`DEVICE`                 |                                                            |:code:`__device__`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL`                  |:code:`__local`                                             |:code:`__shared__`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_ARG`              |:code:`__local`                                             |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL`                 |:code:`__global`                                            |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`RESTRICT`               |:code:`restrict`                                            |:code:`__restrict__`                        |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_ID`               |:code:`get_local_id(0)`                                     |:code:`threadIdx.x`                         |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`LOCAL_SIZE`             |:code:`get_local_size(0)`                                   |:code:`blockDim.x`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL_ID`              |:code:`get_global_id(0)`                                    |:code:`(blockIdx.x*blockDim.x+threadIdx.x)` |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GLOBAL_SIZE`            |:code:`get_global_size(0)`                                  |:code:`(blockDim.x*gridDim.x)`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`GROUP_ID`               |:code:`get_group_id(0)`                                     |:code:`blockIdx.x`                          |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`NUM_GROUPS`             |:code:`get_num_groups(0)`                                   |:code:`gridDim.x`                           |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`SYNC_THREADS`           |:code:`barrier(CLK_LOCAL_MEM_FENCE+CLK_GLOBAL_MEM_FENCE);`  |:code:`__syncthreads();`                    |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`SYNC_WARPS`             | | if SIMT width >= 32:                                     | | if compute capability >= 7.0:            |
|                               | | :code:`mem_fence(CLK_LOCAL_MEM_FENCE)`                   | | :code:`__syncwarp();`                    |
|                               | | otherwise:                                               | | otherwise empty                          |
|                               | | :code:`barrier(CLK_LOCAL_MEM_FENCE)`                     |                                            |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`MEM_FENCE`              |:code:`mem_fence(CLK_LOCAL_MEM_FENCE+CLK_GLOBAL_MEM_FENCE);`|:code:`__threadfence_block();`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+
|:code:`ATOMIC_ADD(dest, value)`|:code:`atom_add(dest, value)`                               |:code:`atomicAdd(dest, value)`              |
+-------------------------------+------------------------------------------------------------+--------------------------------------------+

A few other symbols may or may not be defined based on the device you are running on:
:code:`SUPPORTS_DOUBLE_PRECISION` and :code:`SUPPORTS_64_BIT_ATOMICS`\ .  You
can use :code:`#ifdef` blocks with these symbols to conditionally compile code
based on the features supported by the device.  In addition, the CUDA compiler
defines the symbol :code:`__CUDA_ARCH__`\ , so you can check for this symbol if
you want to have different code blocks for CUDA and OpenCL.

Both OpenCL and CUDA define vector types like :code:`int2` and :code:`float4`\ .
The types they support are different but overlapping.  When writing common code,
use only the vector types that are supported by both OpenCL and CUDA: 2, 3, and 4
element vectors of type :code:`short`\ , :code:`int`\ , :code:`float`\ , and
:code:`double`\ .

CUDA uses functions to construct vector values, such as :code:`make_float2(x, y)`\ .
OpenCL instead uses a typecast like syntax: :code:`(float2) (x, y)`\ .  In common
code, use the CUDA style :code:`make_` functions.  OpenMM provides definitions
of these functions when compiling as OpenCL.

In CUDA, vector types are simply data structures.  You can access their elements,
but not do much more with them.  In contrast, OpenCL's vectors are mathematical
types.  All standard math operators are defined for them, as well as geometrical
functions like :code:`dot()` and :code:`cross()`\ .  When compiling kernels as
CUDA, OpenMM provides definitions of these operators and functions.

OpenCL also supports "swizzle" notation for vectors.  For example, if :code:`f`
is a :code:`float4` you can construct a vector of its first three elements
by writing :code:`f.xyz`\ , or you can swap its first two elements by writing
:code:`f.xy = f.yx`\ .  Unfortunately, there is no practical way to support this
in CUDA, so swizzle notation cannot be used in common code.  Because stripping
the final element from a four component vector is such a common operation, OpenMM
provides a special function for doing it: :code:`trimTo3(f)` is a vector of its
first three elements.

64 bit integers are another data type that needs special handling.  Both OpenCL
and CUDA support them, but they use different names for them: :code:`long` in OpenCL,
:code:`long long` in CUDA.  To work around this inconsistency, OpenMM provides
the typedefs :code:`mm_long` and :code:`mm_ulong` for signed and unsigned 64 bit
integers in device code.

Writing Host Code
*****************

Host code for Common Compute is very similar to host code for OpenCL or CUDA.
In fact, most of the classes provided by the OpenCL and CUDA platforms are
subclasses of Common Compute classes.  For example, OpenCLContext and
CudaContext are both subclasses of ComputeContext.  When writing common code,
each KernelImpl should expect a ComputeContext to be passed to its constructor.
By using the common API provided by that abstract class, it can be used for
either OpenCL or CUDA just based on the particular context passed to it at
runtime.  Similarly, OpenCLNonbondedUtilities and CudaNonbondedUtilities are
subclasses of the abstract NonbondedUtilities class, and so on.

ArrayInterface is an abstract class defining the interface for arrays stored on
the device.  OpenCLArray and CudaArray are both subclasses of it.  To simplify
code that creates and uses arrays, there is also a third subclass called
ComputeArray.  It acts as a wrapper around an OpenCLArray or CudaArray,
automatically creating an array of the appropriate type for the current
platform.  In practice, just follow these rules:

  1. Whenever you need to create an array, make it a ComputeArray.

  2. Whenever you write a function that expects an array to be passed to it,
     declare the type to be ArrayInterface.

If you do these two things, all differences between platforms will be handled
automatically.

OpenCL and CUDA have quite different APIs for compiling and invoking kernels.
To hide these differences, OpenMM provides a set of abstract classes.  To compile
device code, pass the source code to :code:`compileProgram()` on the ComputeContext.
This returns a ComputeProgram.  You can then call its :code:`createKernel()`
method to get a ComputeKernel object, which has methods for setting arguments
and invoking the kernel.

Sometimes you need to refer to vector types in host code, such as to set the
value for a kernel argument or to access the elements of an array.  OpenCL and
CUDA both define types for them, but they have different names, and in any case
you want to avoid using OpenCL-specific or CUDA-specific types in common code.
OpenMM therefore defines types for vectors in host code.  They have the same
names as the corresponding types in device code, only with the prefix :code:`mm_`\ ,
for example :code:`mm_int2` and :code:`mm_float4`\ .

Three component vectors need special care in this context, because the platforms
define them differently.  In OpenCL, a three component vector is essentially a
four component vector whose last component is ignored.  For example,
:code:`sizeof(float3)` is 12 in CUDA but 16 in OpenCL.  Within a kernel this
distinction can usually be ignored, but when communicating between host and
device it becomes vitally important.  It is generally best to avoid storing
three component vectors in arrays or passing them as arguments.  There are no
:code:`mm_` host types defined for three component vectors, because CUDA and
OpenCL would require them to be defined in different ways.
