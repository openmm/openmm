export PATH=$HOME/miniconda/bin:$PATH
export SWIG_LIB=$HOME/miniconda/share/swig/
export CC="clang++"

git clone -b vagrant https://github.com/simtk/openmm.git
cd openmm
conda install --file tools/ci/requirements-conda.txt --yes
conda build tools/conda-recipe
