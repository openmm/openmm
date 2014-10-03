export PATH=$HOME/miniconda/bin:$PATH
export SWIG_LIB=$HOME/miniconda/share/swig/
export CC="clang++"

git clone https://github.com/simtk/openmm.git
cd openmm
git checkout tags/6.1  # To checkout specific release for packaging. 
conda install --file devtools/ci/requirements-conda.txt --yes

# For CI build:
conda build devtools/conda-recipe

# For release build:
cd ../
git clone https://github.com/omnia-md/conda-recipes.git
conda build conda-recipes/openmm

# To upload the file, do something the following command but with the package version changed:

binstar upload -u omnia /home/vagrant/miniconda/conda-bld/linux-64/openmm-6.1-py27_0.tar.bz2
