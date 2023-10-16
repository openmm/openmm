.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++


.. _the-cpu-platform:

The CPU Platform
################

CpuPlatform is a subclass of ReferencePlatform.  It provides optimized versions
of a small number of kernels, while using the reference implementations for all
the others.  Any kernel implementation written for the reference Platform will
work equally well with the CPU platform.  Of course, if that kernel happens to
be a performance bottleneck, you will probably want to write an optimized
version of it.  But many kernels have negligible effect on performance, and for
these you can just use the same implementation for both platforms.

If you choose to do that, you can easily support both platforms with a single
plugin library.  Just implement :code:`registerKernelFactories()` like this:
::

    extern "C" void registerKernelFactories() {
        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
            Platform& platform = Platform::getPlatform(i);
            if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
                // Create and register your KernelFactory.
            }
        }
    }

The loop identifies every ReferencePlatform, either an instance of the base
class or of a subclass, and registers a KernelFactory for every one.
