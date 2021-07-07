## How to Contribute to OpenMM Development

We welcome anyone who wants to contribute to the project, whether by adding a feature,
fixing a bug, or improving documentation.  The process is quite simple.

First, it is always best to begin by opening an issue on GitHub that describes the change you
want to make.  This gives everyone a chance to discuss it before you put in a lot of work.
For bug fixes, we will confirm that the behavior is actually a bug and that the proposed fix
is correct.  For new features, we will decide whether the proposed feature is something we
want and discuss possible designs for it.

Once everyone is in agreement, the next step is to
[create a pull request](https://help.github.com/en/articles/about-pull-requests) with the code changes.
For larger features, feel free to create the pull request even before the implementation is
finished so as to get early feedback on the code.  When doing this, put the letters "WIP" at
the start of the title of the pull request to indicate it is still a work in progress.

For new features, consult the [New Feature Checklist](https://github.com/openmm/openmm/wiki/Checklist-when-adding-a-new-feature),
which lists various items that need to be included before the feature can be merged (documentation,
tests, serialization, support for all APIs, etc.).  Not every item is necessarily applicable to
every new feature, but usually at least some of them are.

The core developers will review the pull request and may suggest changes.  Simply push the
changes to the branch that is being pulled from, and they will automatically be added to the
pull request.  In addition, the full test suite is automatically run on every pull request,
and rerun every time a change is added.  Once the tests are passing and everyone is satisfied
with the code, the pull request will be merged.  Congratulations on a successful contribution!

## Building OpenMM

OpenMM uses [cmake](https://cmake.org/) as our build system and requires quite a number of dependencies. Detailed instructions for compiling OpenMM from source are found in the [User Guide](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code). The following is a summary of the most important steps to compile the CPU platform and documentation.

### Installing dependencies with Conda

We provide the core dependencies for the CPU platform and documentation in a Conda environment file. They can be installed easily:

```shell
conda env create -f devtools/environment.yml
```

Since this environment includes the compiler toolchains necessary to build OpenMM, the build system will not be able to use any system-installed GPU SDK. On NVIDIA GPUs, the CUDA toolkit can be installed relatively easily:

```shell
conda install -n openmm-dev -c conda-forge cudatoolkit-dev
```

Unfortunately, no equivalent metapackage exists for OpenCL.
On 64 bit x86 platforms, try the `mesa-libgl-devel-cos6-x86_64` package; otherwise, check the [available Mesa development packages on `conda-forge`](https://anaconda.org/search?q=access%3Apublic+type%3Aconda+mesa-libgl-devel) for your architecture.
You may also need an OCL loader; try `ocl-icd-system` on Linux or `khronos-opencl-icd-loader` on Windows or MacOS.

### Configuring make files

```shell
mkdir build
cd build
# Activate the Conda environment
conda activate openmm-dev
# Configure build. For other options, check the User Guide or try `ccmake ..`
# Remove the last switch to try building the OpenCL platform
cmake .. \
    -DOPENMM_GENERATE_API_DOCS=ON \
    -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
    -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DOPENMM_BUILD_PYTHON_WRAPPERS=ON \
    -DOPENMM_BUILD_OPENCL_LIB=OFF
```

### Compiling OpenMM

From the `build` directory created in the previous step:

```shell
# Activate the Conda environment
conda activate openmm-dev
# Build libraries and tests
make
# Build the Python layer and install the OpenMM libraries to the conda environment
make install PythonInstall
```

### Compiling the documentation

```shell
# Activate the Conda environment
conda activate openmm-dev
# Compile the User's and Developer's guides as HTML
make sphinxhtml
# Compile the C++ API docs
make C++ApiDocs
# Compile the Python API docs
make PythonApiDocs
```
