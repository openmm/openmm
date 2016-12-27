#!/bin/bash
set -e -x

pushd .
cd ..

# Update homebrew
brew uninstall -y brew-cask || brew untap -y caskroom/cask || 1
brew update -y --quiet
brew tap -y caskroom/cask

# Install doxygen
brew install -y --quiet doxygen

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
sudo tlmgr install titlesec framed threeparttable wrapfig multirow collection-fontsrecommended hyphenat xstring fncychap tabulary
ls -ltr /usr/local/texlive/*/texmf-dist/tex/latex/

# Build packages
source openmm/devtools/packaging/scripts/osx/prepare.sh
source openmm/devtools/packaging/scripts/osx/build.sh
source openmm/devtools/packaging/scripts/osx/package.sh

# Copy packages
cp -r packaging /openmm/

popd
