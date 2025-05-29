.. role:: code
.. raw:: html

    <style> .code {font-family:monospace;} </style>
    <style> .caption {text-align:center;} </style>

.. highlight:: c++

.. _writing-plugins:

Writing Plugins
###############

A plugin is a dynamic library that adds new features to OpenMM.  It is typically
stored in the :code:`lib/plugins` directory inside your OpenMM installation,
and gets loaded along with all other plugins when the user calls
::

    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

It is also possible to load plugins from a different directory, or to load them
individually by calling :code:`Platform::loadPluginLibrary()`\ .

Every plugin must implement two functions that are declared in the
PluginInitializer.h header file:
::

    extern "C" void registerPlatforms();
    extern "C" void registerKernelFactories();

When a plugin is loaded, these two functions are invoked to register any
Platforms and KernelFactories defined by the plugin.  When many plugins are
loaded at once by calling :code:`Platform::loadPluginsFromDirectory()`\ ,
:code:`registerPlatforms()` is first called on all of them, then
:code:`registerKernelFactories()` is called on all of them.  This allows one
plugin to define a Platform, and a different plugin to add KernelFactories to
it; the Platform is guaranteed to be registered by the first plugin before the
second plugin tries to add its KernelFactories, regardless of what order the
plugins happen to be loaded in.

Creating New Platforms
**********************

One common type of plugin defines a new Platform.  There are four such plugins
that come with OpenMM: one for the Reference platform, one for the CPU Platform,
one for the CUDA Platform, and one for the OpenCL Platform.

To define a new Platform, you must create subclasses of the various abstract
classes in the OpenMM Low Level API: a subclass of Platform, one or more
subclasses of KernelFactory, and a subclass of each KernelImpl.  That is easy to
say, but a huge amount of work to actually do.  There are many different
algorithms involved in computing forces, enforcing constraints, performing
integration, and so on, all of which together make up a Platform.  Of course,
there is no requirement that every Platform must implement every possible
feature.  If you do not provide an implementation of a particular kernel, it
simply means your Platform cannot be used for any simulation that requires that
kernel; if a user tries to do so, an exception will be thrown.

Your plugin’s :code:`registerPlatforms()` function should create an instance
of your Platform subclass, then register it by calling
:code:`Platform::registerPlatform()`\ .  You also must register the
KernelFactory for each kernel your Platform supports.  This can be done in the
:code:`registerKernelFactories()` function, or more simply, directly in the
Platform’s constructor.  You can use as many different KernelFactories as you
want for different kernels, but usually it is simplest to use a single
KernelFactory for all of them.  The support for multiple KernelFactories exists
primarily to let plugins add new features to existing Platforms, as described in
the next section.

Creating New Forces
*******************

Another common type of plugin defines new Forces and provides implementations of
them for existing Platforms.  (Defining new Integrators is not specifically
discussed here, but the process is very similar.)  There are two such plugins
that come with OpenMM.  They implement the AMOEBA force field and Drude
oscillators, respectively.

As an example, suppose you want to create a new Force subclass called
StringForce that uses the equations of String Theory to compute the interactions
between particles.  You want to provide implementations of it for all four
standard platforms: Reference, CPU, CUDA, and OpenCL.

The first thing to realize is that this *cannot* be done with only a plugin
library.  Plugins are loaded dynamically at runtime, and they relate to the low
level API; but you must also provide a public API.  Users of your class need to
create StringForce objects and call methods on them.  That means providing a
header file with the class declaration, and a (non-plugin) library with the
class definition to link their code against.  The implementations for particular
Platforms can be in plugins, but the public API class itself cannot.  Or to put
it differently, the full “plugin” (from the user’s perspective) consists of
three parts: the library OpenMM loads at runtime (which is what OpenMM considers
to be the “plugin”), a second library for users to link their code against, and
a header file for them to include in their source code.

To define the API, you will need to create the following classes:

#. StringForce.  This is the public API for your force, and users will directly
   link against the library containing it.
#. StringForceImpl.  This is the ForceImpl subclass corresponding to
   StringForce.  It should be defined in the same library as StringForce, and
   StringForce’s :code:`createImpl()` method should create an instance of it.
#. CalcStringForceKernel.  This is an abstract class that extends KernelImpl,
   and defines the API by which StringForceImpl invokes its kernel.  You only need
   to provide a header file for it, not an implementation; those will be provided
   by Platforms.


Now suppose you are writing the OpenCL implementation of StringForce.  Here are
the classes you need to write:

#. OpenCLCalcStringForceKernel.  This extends CalcStringForceKernel and provides
   implementations of its virtual methods.  The code for this class will probably
   be very complicated (and if it actually works, worth a Nobel Prize).  It may
   execute many different GPU kernels and create its own internal data structures.
   But those details are entirely internal to your own code.  As long as this class
   implements the virtual methods of CalcStringForceKernel, you can do anything you
   want inside it.
#. OpenCLStringForceKernelFactory.  This is a KernelFactory subclass that knows
   how to create instances of OpenCLCalcStringForceKernel.


Both of these classes should be packaged into a dynamic library (.so on Linux,
.dylib on Mac, .dll on Windows) that can be loaded as a plugin.  This library
must also implement the two functions from PluginInitializer.h.
:code:`registerPlatforms()` will do nothing, since this plugin does not
implement any new Platforms.  :code:`registerKernelFactories()` should call
\ :code:`Platform::getPlatform("OpenCL")` to get the OpenCL Platform,
then create a new OpenCLStringForceKernelFactory and call
:code:`registerKernelFactory()` on the Platform to register it.  If the OpenCL
Platform is not available, you should catch the exception then return without
doing anything.  Most likely this means there is no OpenCL runtime on the
computer your code is running on.
