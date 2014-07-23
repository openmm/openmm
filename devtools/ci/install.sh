# Original OpenMM .travis.yml requirements.
sudo apt-get update -qq
sudo apt-get install -qq libpcre3 libpcre3-dev gromacs
sudo apt-get install -qq swig doxygen llvm-3.3
sudo apt-get install -qq libgl1-mesa-dev opencl-headers fglrx=2:8.960-0ubuntu1 # for opencl support
export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.3

# New requirements.
sudo apt-get install -qq -y g++ gfortran csh
sudo apt-get install -qq -y g++-multilib gcc-multilib
wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh
bash Miniconda-3.0.5-Linux-x86_64.sh -b
PIP_ARGS="-U"

export PATH=$HOME/miniconda/bin:$PATH

conda update --yes conda
conda config --add channels http://conda.binstar.org/omnia
conda create --yes -n ${python} python=${python} --file tools/ci/requirements-conda.txt
source activate $python
$HOME/miniconda/envs/${python}/bin/pip install $PIP_ARGS nose-exclude
