<!-- Authored by Jaime Rodríguez-Guerra, Chodera Lab. December 2020 -->

# Our Continuous Integration setup

OpenMM can be described as a C++ library with wrappers available in different programming languages (Python, C, Fortran). The heavy lifting is performed by the backend platforms, which can be based on CPU, CUDA and/or OpenCL (and possibly more in the future). All of this is supported for different operating systems and architectures. As a result, the CI setup can get a bit involved, but this document will try to clarify how it works and what we support.

## Implementation overview

OpenMM's CI runs mainly on GitHub Actions, with one separate Jenkins box running the GPU tests (generously provided by Jason Swails).

The build matrix covers:

- Operating systems and architecture:
  - Linux x64
  - MacOS Intel
  - Windows
  - Linux ppc64le (PowerPC)
  - Linux aarch64 (ARM)
- Python
  - CPython 3.6, 3.7, 3.8, 3.9
- CUDA versions
  - 10.0 and above (Linux x64, Linux ppc64le, Windows)
- OpenCL implementations
  - Nvidia (tested along CUDA)
  - AMD 3.0
- Sysroots and C++ Compilers
  - Linux: System's GCC 7 and whatever conda-forge is pinning (GCC 9 as of writing)
  - MacOS: System's, targetting 10.9 SDK
  - Windows: VS2019

Before I describe the pipelines, I will clarify some concepts and idiosyncrasies in GitHub Actions

- The configuration file lives on `.github/workflows/CI.yml`. This directory can host more than one YML _workflow_, each describing a set of event that will trigger a run.
- The workflow specifies a set of triggers (key `on`) and a list of `jobs` to run. We run the `CI` workflow for:
  - Pushes to `master`
  - Pull requests targetting `master`
  - Nightlies
- Currently, the workflow contains four jobs: `unix`, `windows`, `docker`, `docs`. Each job can be run several times, depending on the configuration of `jobs.*.strategy.matrix`. All those jobs replicas will run in parallel and individually. The [`Actions > Summary`](https://github.com/openmm/openmm/actions/runs/451301350) overview can help visualize this.
- Within each job, you find `steps`. A step can either run a script on a `shell` or use a GitHub `action` to perform a task.
  - For example, cloning the repo or setting up Miniconda are both independent GitHub _actions_. You will recognize this because they contain the keyword `uses:`.
  - Running CMake is a shell step, which uses `run:`.
  - Note 1: Each step is run a new shell session. Environment variables won't survive across steps, unless you add them to the `$GITHUB_ENV` file: `echo "VARIABLE=VALUE" >> ${GITHUB_ENV}`. You can also use step `outputs` but that's more involved and rarely needed.
  - Note 2: Due to the design of `conda-incubator/setup-miniconda`, all subsequent steps that rely on a conda environment require us to specify an OS-dependent custom shell. Do remember this if you need to add more steps in the job!
- Steps can be run or skipped based on conditions expressed inside an `if:` key. This is how we control whether we need to install CUDA or not, for example. Jobs can have `if` check, if needed.
- Steps can define environment variables in their `env:` key, but they will only be available in that step. A `job` can do it too, and these will be available for all steps.

## Details per operating system

The different implementations are very similar to what we do on Linux x64, so I will explain this one on detail and the rest will only comment on the relevant differences.

### Linux x64

- Part of the `unix` pipeline.
- Runs on `ubuntu-latest`, as provided by GitHub Actions.
- Uses `conda-incubator/setup-miniconda` to setup the bundled Miniconda and install a conda environment available that provides the building and testing dependencies (CMake, Swig, the adequate Python version, etc). These environment files are located under `devtools/ci/gh-actions/conda-envs`, per operating system.
- Depending on the matrix configuration, we also install CUDA and/or AMD's OpenCL. These conditional steps are evaluated using GHA's builtin `if` mechanism. Ideally we would install this within the conda environment, but sometimes they are not available (licensing issues, etc(), so we delegate that to the system packages or vendor installers.
  - For CUDA, we check whether `cuda-version` is not empty, and pass it to `devtools/ci/gh-actions/scripts/install_cuda.sh` as an environment variable.
  - For OpenCL, we check whether `OPENCL` is `true` and run `devtools/ci/gh-actions/scripts/install_amd_opencl.sh`. This relies on a installer located in a S3 bucket. This could be refactored to install different OpenCL implementations (ROCm, Intel, etc).
- Some matrix entries require us to install the conda forge compilers, which are used instead of the system's if present.
- Now we need to configure and download the CCache contents. The keys are built off the matrix name, and a `YYYYDDMM-HHMMSS` timestamp. A secret `CACHE_VERSION` is also included so one can bump the cache by modifying this secret in the repository settings. The configuration is done through environment variables defined at the beginning of the job (key `jobs.unix.env`).
- CMake is finally invoked, targetting the conda environment as destination (`CONDA_PREFIX`). Additional flags are passed from the matrix configuration. This is how we enable or disable features per matrix entry.
- CCache performance is assessed.
- Then we build the C++ libraries and Python wrappers, but separately. This way we can visually check which part failed more easily. Tests are also run separately for the same reason. Whether Python is built and/or tested is checked through the contents of `CMAKE_FLAGS`.

### MacOS Intel

- Part of the `unix` pipeline.
- Runs on `macos-latest`.
- Uses `conda-incubator/setup-miniconda`, pointing to the relevant environment file.
- Neither CUDA nor OpenCL installation scripts are run. Instead, we download and install the 10.9 SDK using `devtools/ci/gh-actions/scripts/install_macos_sdk.sh`. This is done so we can mimic what Conda Forge does in their feedstocks. Check the scripts comments for more info.
- Everything else is the same.

### Windows

- Sole member of the `windows` pipeline.
- Runs on `windows-latest`.
- Uses `conda-incubator/setup-miniconda`, pointing to the relevant environment file.
- Installs CUDA with the Nvidia installers using `devtools/ci/gh-actions/scripts/install_cuda.bat`, which requires an environment variable `CUDA_VERSION`, exported from the corresponding matrix entry. Again, this only runs if `matrix.cuda-version` is not empty.
- Everything else is the same.

### PowerPC & ARM

- Part of the `docker` pipeline.
- These run on a Docker image on top of `ubuntu-latest`. The Docker image itself depends on the architecture chosen (ppc64le, aarch64) and what CUDA version we want. These are provided by Conda Forge, so they have `conda` preinstalled and ready to go.
- Since it's a different architecture, we need to configure QEMU first. This is done automatically with a Docker image, mimicking what Conda Forge does.
- We start the Docker image. The working directory (`$GITHUB_WORKSPACE`) is mounted with read/write permissions on `/home/conda/workspace`, so we can communicate back with the host using files, and also use CCache.
- The Docker image will run `devtools/ci/gh-actions/scripts/run_steps_inside_docker_image.sh`. This script mostly does what you saw for Linux x64, with some differences:
  - We don't need to install CUDA or setup Miniconda, because they are preinstalled in the Docker image.
  - We patch some dependencies from the environment file because they are not available for this architecture. To save one conda environment solve, we also patch the Python version in the environment file.
  - These images don't come with a system compiler, so we specify one in the matrix configuration:
    - If `compilers` contains a value that starts with `devtoolset-`, we understand we want a CentOS devtoolse. So far, we specify `devtoolset-7`.
    - If `compilers` is any other thing, we understand that's a (space-separated series of) Conda packages. Since Conda Forge provides a metapackage named `compilers` that will install all of them for the current platform, we use that one. That's why some entries have a `compilers: compilers` entry.
  - Everything else runs as usual.
- Do note that the whole Docker run is a single GitHub Actions step, so it's not as visually appealing. I tried my best to group the commands with the `::group::` syntax so it's easier to follow, but it's not the same.
- If the script runs successfully, it will create an empty file. We test for existence after the Docker run to make sure.

> Note: Since these use software emulation, they are really slow. Still, they can run successfully within the 6h GHA provides. If GHA upgrades to better CI machines with hardware based virtualization, they might be able to run with close-to-native performance.

### Docs

This is a Linux-x64 pipeline optimized for building the documentation only. It's provided as a separate entry because I didn't want to overcomplicate the `if:` logic in the `unix` pipeline. It's essentially the same, but:

- It uses a different environment file in `setup-miniconda`.
- It only builds the docs, and their dependencies. No tests, for example.
- It contains a deployment step, which will copy the contents to the S3 bucket _only_ when run on `master`, ignoring cron jobs. The required secrets must be defined in the repository settings with the following exact key names. Just copy paste the values there. GitHub will encrypt and mask them.
  - `AWS_S3_BUCKET`
  - `AWS_ACCESS_KEY_ID`
  - `AWS_SECRET_ACCESS_KEY`
- It will also check for dead links using a Node package. This is run _after_ deployment so it won't prevent that, but it will still signal the job as failed if the docs contain broken links.

## Shortcomings

There are some limitations when compared to other CI services, but I guess this list will be shorter over time:

- Cache cannot be invalidated directly. Instead, I included a secret `CACHE_VERSION` that is part of the cache key. If you change the value of this secret, it will functionally prevent access to the previous cache. It also expires every 7 days. Note that since this trick uses a secret, the value of `CACHE_VERSION` will be masked in the log output. As a result, make sure to use something short but meaningless and difficult to find in the wild (e.g. `pqgbhl` instead of `0`).
- There's no `ci skip` functionality (yet).

## Extra content

### How to debug PowerPC / ARM locally

From the root of the repository, run the following script. There are
some variables you might want to edit (PPC vs ARM, Python version, etc).
Take a look to the script first in that case.

```bash
bash devtools/ci/gh-actions/start_docker_locally.sh
```

You will be inside the Docker image after a few moments. The repo root has
been mounted to `/home/conda/workspace`.

Run this other script to reproduce the CI steps exactly. Do NOT `source` scripts,
since a failure will exit Docker altogether. Always use new `bash` processes
to avoid starting from scratch.

```bash
bash /home/conda/workspace/devtools/ci/gh-actions/scripts/run_steps_inside_docker_image.sh
```
