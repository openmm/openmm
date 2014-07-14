# Prepare a vagrant CentOS 6.5 VM for building OpenMM
# Needs latest version of vagrant to auto-download the chef package
#vagrant init chef/centos-6.5
#vagrant up
#vagrant ssh


# Download and enable the EPEL RedHat EL extras repository
mkdir ~/Software
cd Software
sudo yum install wget -y
wget http://mirror.umd.edu/fedora/epel/6/i386/epel-release-6-8.noarch.rpm
sudo rpm -i epel-release-6-8.noarch.rpm

sudo yum update -y

# Several of these come from the EPEL repo
sudo yum install clang-3.4 cmake28 graphviz perl flex bison rpm-build texlive texlive-latex ghostscript gcc gcc-c++ git vim -y

# Probably can't use RHEL6 version of doxygen because it's very old.
wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.7.src.tar.gz
rpmbuild -ta doxygen-1.8.7.src.tar.gz
sudo rpm -i ~/rpmbuild/RPMS/x86_64/doxygen-1.8.7-1.x86_64.rpm
rm ~/rpmbuild -r


sudo yum clean headers
sudo yum clean packages

# Install CUDA6 for RHEL6
cd ~/Software
wget http://developer.download.nvidia.com/compute/cuda/repos/rhel6/x86_64/cuda-repo-rhel6-6.0-37.x86_64.rpm
sudo rpm -i  cuda-repo-rhel6-6.0-37.x86_64.rpm
sudo yum clean expire-cache
sudo yum install cuda -y
rm cuda-repo-rhel6-6.0-37.x86_64.rpm


# Install Conda
cd ~/Software
wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh
bash Miniconda-3.0.5-Linux-x86_64.sh -b

# So there is a bug in some versions of anaconda where the path to swig files is HARDCODED.  Below is workaround.  See https://github.com/ContinuumIO/anaconda-issues/issues/48
sudo ln -s  ~/miniconda/ /opt/anaconda1anaconda2anaconda3

export PATH=$HOME/miniconda/bin:$PATH
conda config --add channels http://conda.binstar.org/omnia
conda install --yes fftw3f jinja2 swig sphinx conda-build cmake


# Download AMD APP SDK from here, requires click agreement: http://developer.amd.com/amd-license-agreement-appsdk/
# Ideally we could cache this on AWS or something...
mkdir ~/Software/AMD
cd ~/Software/AMD
# Copy the tarball to this directory from wherever you got it.
cp /vagrant/AMD-APP-SDK-v2.9-lnx64.tgz  ./
tar -zxvf  /vagrant/AMD-APP-SDK-v2.9-lnx64.tgz
sudo ./Install-AMD-APP.sh

