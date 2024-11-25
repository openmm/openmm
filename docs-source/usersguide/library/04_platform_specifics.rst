.. _platform-specific-properties:

Platform-Specific Properties
############################

When creating a Context, you can specify values for properties specific to a
particular Platform.  This is used to control how calculations are done in ways
that are outside the scope of the generic OpenMM API.

To do this, pass both the Platform object and a map of property values to the
Context constructor:

.. code-block:: c

    Platform& platform = Platform::getPlatform("OpenCL");
    map<string, string> properties;
    properties["DeviceIndex"] = "1";
    Context context(system, integrator, platform, properties);

After a Context is created, you can use the Platform’s \
:code:`getPropertyValue()` method to query the values of properties.

OpenCL Platform
***************

The OpenCL Platform recognizes the following Platform-specific properties:

* Precision: This selects what numeric precision to use for calculations.
  The allowed values are “single”, “mixed”, and “double”.  If it is set to
  “single”, nearly all calculations are done in single precision.  This is the
  fastest option but also the least accurate.  If it is set to “mixed”, forces are
  computed in single precision but integration is done in double precision.  This
  gives much better energy conservation with only a slight decrease in speed.
  If it is set to “double”, all calculations are done in double precision.  This
  is the most accurate option, but is usually much slower than the others.
* UseCpuPme: This selects whether to use the CPU-based PME
  implementation.  The allowed values are “true” or “false”.  Depending on your
  hardware, this might (or might not) improve performance.
* OpenCLPlatformIndex: When multiple OpenCL implementations are installed on
  your computer, this is used to select which one to use.  The value is the
  zero-based index of the platform (in the OpenCL sense, not the OpenMM sense) to use,
  in the order they are returned by the OpenCL platform API.  This is useful, for
  example, in selecting whether to use a GPU or CPU based OpenCL implementation.
* DeviceIndex: When multiple OpenCL devices are available on your
  computer, this is used to select which one to use.  The value is the zero-based
  index of the device to use, in the order they are returned by the OpenCL device
  API.


The OpenCL Platform also supports parallelizing a simulation across multiple
GPUs.  To do that, set the DeviceIndex property to a comma separated list
of values.  For example,

.. code-block:: c

    properties["DeviceIndex"] = "0,1";

This tells it to use both devices 0 and 1, splitting the work between them.

CUDA Platform
*************

The CUDA Platform recognizes the following Platform-specific properties:

* Precision: This selects what numeric precision to use for calculations.
  The allowed values are “single”, “mixed”, and “double”.  If it is set to
  “single”, nearly all calculations are done in single precision.  This is the
  fastest option but also the least accurate.  If it is set to “mixed”, forces are
  computed in single precision but integration is done in double precision.  This
  gives much better energy conservation with only a slight decrease in speed.
  If it is set to “double”, all calculations are done in double precision.  This
  is the most accurate option, but is usually much slower than the others.
* UseCpuPme: This selects whether to use the CPU-based PME implementation.
  The allowed values are “true” or “false”.  Depending on your hardware, this
  might (or might not) improve performance.
* TempDirectory: This specifies a directory where temporary files can be
  written while compiling kernels.  OpenMM usually can locate your operating
  system’s temp directory automatically (for example, by looking for the TEMP
  environment variable), so you rarely need to specify this.
* DeviceIndex: When multiple CUDA devices are available on your computer,
  this is used to select which one to use.  The value is the zero-based index of
  the device to use, in the order they are returned by the CUDA API.
* UseBlockingSync: This is used to control how the CUDA runtime
  synchronizes between the CPU and GPU.  If this is set to “true” (the default),
  CUDA will allow the calling thread to sleep while the GPU is performing a
  computation, allowing the CPU to do other work.  If it is set to “false”, CUDA
  will spin-lock while the GPU is working.  Setting it to "false" can improve performance slightly,
  but also prevents the CPU from doing anything else while the GPU is working.
* DeterministicForces: In some cases, the CUDA platform may compute forces
  in ways that are not fully deterministic (typically differing in what order a
  set of numbers get added together).  This means that if you compute the forces
  twice for the same particle positions, there may be tiny differences in the
  results.  In most cases this is not a problem, but certain algorithms depend
  on forces being exactly reproducible to the last bit.  If you set this
  property to "true", it will instead do these calculations in a way that
  produces fully deterministic results, at the cost of a small decrease in
  performance.

The CUDA Platform also supports parallelizing a simulation across multiple GPUs.
To do that, set the DeviceIndex property to a comma separated list of
values.  For example,

.. code-block:: c

    properties["DeviceIndex"] = "0,1";

This tells it to use both devices 0 and 1, splitting the work between them.

HIP Platform
************

The HIP Platform recognizes exactly the same Platform-specific properties as
the CUDA platform.

CPU Platform
************

The CPU Platform recognizes the following Platform-specific properties:

* Threads: This specifies the number of CPU threads to use.  If you do not
  specify this, OpenMM will select a default number of threads as follows:

  * If an environment variable called OPENMM_CPU_THREADS is set, its value is
    used as the number of threads.
  * Otherwise, the number of threads is set to the number of logical CPU cores
    in the computer it is running on.

  Usually the default value works well.  This is mainly useful when you are
  running something else on the computer at the same time, and you want to
  prevent OpenMM from monopolizing all available cores.

.. _platform-specific-properties-determinism:

Determinism
***********

Whether a simulation is deterministic will depend on what platform you run on in
addition to what settings/methods you use. For instance, as of this writing,
using PME on the Reference, OpenCL, and double-precision CUDA will result in
deterministic simulations. Single-precision CUDA and CPU platforms are not
deterministic in this case. However, none of this behavior is guaranteed in
future versions. In many cases it will still result in an identical trajectory.
If determinism is a critical for your needs, you should carefully check to
ensure that your settings and platform allow for this.
