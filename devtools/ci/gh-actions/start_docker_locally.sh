#!/bin/bash

# This is an example script on how to debug locally with Docker!
# If it does not work, it might be out of date. In that case,
# check the steps used in /.github/workflows/CI.yml

set -euxo pipefail

# This is the image for PowerPC + CUDA
export DOCKER_IMAGE="quay.io/condaforge/linux-anvil-ppc64le-cuda:10.2"
# # Use this other one for ARM debugging
# export DOCKER_IMAGE="quay.io/condaforge/linux-anvil-aarch64"

# With Conda Forge compilers (GCC9)
export COMPILERS="compilers"
# # With RH devtoolset (GCC7)
# export COMPILERS="devtoolset-7"

# Choose your Python version
export PYTHON_VER="3.9"

# Number of CPUs to use
export CPU_COUNT=2

echo "Preparing Docker..."
docker run --rm --privileged multiarch/qemu-user-static:register --reset --credential yes
ls /proc/sys/fs/binfmt_misc/

docker info

# In order for the conda-build process in the container to write to the mounted
# volumes, we need to run with the same id as the host machine, which is
# normally the owner of the mounted volumes, or at least has write permission
export HOST_USER_ID=$(id -u)

# Check if docker-machine is being used (normally on OSX) and get the uid from
# the VM
if hash docker-machine 2> /dev/null && docker-machine active > /dev/null; then
    export HOST_USER_ID=$(docker-machine ssh $(docker-machine active) id -u)
fi

docker run \
  -it \
  -v "$(pwd)":/home/conda/workspace:rw,z \
  -e HOST_USER_ID \
  -e CPU_COUNT \
  -e PYTHON_VER \
  -e COMPILERS \
  ${DOCKER_IMAGE} \
  bash

# Once you are inside the Docker session, you can use this to reproduce the CI steps:
#
# bash /home/conda/workspace/devtools/ci/gh-actions/scripts/run_steps_inside_docker_image.sh