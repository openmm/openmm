#!/bin/bash

bash -e devtools/ci/jenkins/install.sh \
        -DOPENMM_BUILD_CUDA_LIB=false \
        -DOPENMM_BUILD_OPENCL_LIB=false && \
bash -e devtools/ci/jenkins/test.sh -R 'Test(Cpu|Reference)' --parallel 4
bash -e devtools/ci/jenkins/test_python.sh
