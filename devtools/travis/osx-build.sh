#!/bin/bash
set -e -x

pushd .
cd ..

# Update homebrew
brew uninstall -y brew-cask || brew untap -y caskroom/cask || 1
brew update -y --quiet
brew tap -y caskroom/cask

# Install Miniconda
MINICONDA=Miniconda2-latest-Linux-x86_64.sh
MINICONDA_HOME=$HOME/miniconda
MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget -q https://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b -p $MINICONDA_HOME

# Configure miniconda
export PIP_ARGS="-U"
export PATH=$MINICONDA_HOME/bin:$PATH
conda update --yes conda
conda install --yes conda-build jinja2 anaconda-client pip

# Install doxygen
#brew install -y --quiet doxygen
conda install -c conda-forge doxygen

# Install CUDA
curl -O -s http://developer.download.nvidia.com/compute/cuda/$CUDA_VERSION/Prod/network_installers/mac/x86_64/cuda_mac_installer_tk.tar.gz
curl -O -s http://developer.download.nvidia.com/compute/cuda/$CUDA_VERSION/Prod/network_installers/mac/x86_64/cuda_mac_installer_drv.tar.gz
sudo tar -zxf cuda_mac_installer_tk.tar.gz -C /;
sudo tar -zxf cuda_mac_installer_drv.tar.gz -C /;
rm -f cuda_mac_installer_tk.tar.gz cuda_mac_installer_drv.tar.gz

# Install latex.
brew tap -y --quiet Caskroom/cask;
brew cask install -y --quiet basictex
export PATH="/usr/texbin:${PATH}:/usr/bin"
sudo tlmgr update --self
sudo tlmgr install collection-basic collection-fontsrecommended collection-latex
sudo tlmgr install titlesec framed threeparttable wrapfig multirow collection-fontsrecommended hyphenat xstring fncychap tabulary capt-of needspace eqparbox
ls -ltr /usr/local/texlive/*/texmf-dist/tex/latex/

# Build packages
source openmm/devtools/packaging/scripts/osx/prepare.sh
source openmm/devtools/packaging/scripts/osx/build.sh
source openmm/devtools/packaging/scripts/osx/package.sh

# Copy packages
cp -r packaging /openmm/

popd
