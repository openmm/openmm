# This script install CUDA on Ubuntu-based systems
# It uses the Nvidia repos for Ubuntu 22.04
# Future versions might require an updated repo
# It expects a $CUDA_VERSION environment variable set to major.minor (e.g. 10.0)

set -euxo pipefail

# Enable retrying
echo 'APT::Acquire::Retries "5";' | sudo tee /etc/apt/apt.conf.d/80-retries

sudo apt-key del 7fa2af80
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
CUDA_APT=${CUDA_VERSION/./-}
sudo apt-get install -y \
    libgl1-mesa-dev cuda-compiler-${CUDA_APT} \
    cuda-drivers cuda-driver-dev-${CUDA_APT} \
    cuda-cudart-${CUDA_APT} cuda-cudart-dev-${CUDA_APT} \
    libcufft-${CUDA_APT} libcufft-dev-${CUDA_APT} \
    cuda-nvrtc-${CUDA_APT} cuda-nvrtc-dev-${CUDA_APT} \
    cuda-nvprof-${CUDA_APT} cuda-profiler-api-${CUDA_APT}
sudo apt-get clean

export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION}
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}
export PATH=${CUDA_HOME}/bin:${PATH}

echo "CUDA_HOME=${CUDA_HOME}" >> ${GITHUB_ENV}
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> ${GITHUB_ENV}
echo "PATH=${PATH}" >> ${GITHUB_ENV}