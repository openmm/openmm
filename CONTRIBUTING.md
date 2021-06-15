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

OpenMM uses [cmake](https://cmake.org/) as our build system and requires quite a number of dependencies.

### Installing dependencies

OpenMM's dependencies are described in detail in [the User Guide](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code). The easiest and fastest way to install all required dependencies is to use the [Conda package manager](https://conda.io).

#### on Linux

```shell
# Install build requirements
conda env create -n openmm-dev -f devtools/ci/gh-actions/conda-envs/build-ubuntu-latest.yml
# Add documentation requirements to environment - optional
conda env update -n openmm-dev -f devtools/ci/gh-actions/conda-envs/docs.yml
```

#### on OSX

```shell
# Install build requirements
conda env create -n openmm-dev -f devtools/ci/gh-actions/conda-envs/build-macos-latest.yml
# Add documentation requirements to environment - optional
conda env update -n openmm-dev -f devtools/ci/gh-actions/conda-envs/docs.yml
```

#### on Windows

```shell
# Install build requirements
conda env create -n openmm-dev -f devtools/ci/gh-actions/conda-envs/build-windows-latest.yml
# Add documentation requirements to environment - optional
conda env update -n openmm-dev -f devtools/ci/gh-actions/conda-envs/docs.yml
```

### Configuring make files

```shell
mkdir build
cd build
# For other options, check the User Guide or try `ccmake ..`
cmake .. -DOPENMM_GENERATE_API_DOCS=ON
```

### Compiling OpenMM

From the `build` directory created in the previous step:

```shell
# Build libraries and tests
make
# Build Python API
make PythonInstall
```

### Compiling the documentation

```shell
# Compile the User's and Developer's guides as HTML
make sphinxhtml
# Compile the C++ API docs
make C++ApiDocs
# Compile the Python API docs
make PythonApiDocs
```
