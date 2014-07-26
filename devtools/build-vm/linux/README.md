# Utilities for building OpenMM on a virtual machine

Building the VM with [vagrant](http://vagrantup.com):

* Make sure the latest [vagrant 1.6.3](http://vagrantup.com) is installed.
* Download the AMD APP SDK for Linux 64-bit `AMD-APP-SDK-v2.9-lnx64.tgz` from the [AMD APP SDK site](http://developer.amd.com/tools-and-sdks/opencl-zone/opencl-tools-sdks/amd-accelerated-parallel-processing-app-sdk/) and place it in the current (shared) directory.
* Run `vagrant up` to bring up and provision the VM.
* Build the conda package: `vagrant ssh -c "bash /vagrant/build_conda_package_vagrant.sh"` to build the conda package
* Log in to binstar: `binstar login`, giving username and password
* Set `$BINSTAR_TOKEN` to the appropriate binstar access token: `export BINSTAR_TOKEN=$(binstar auth create --org omnia --name vagrant-test)`
* Upload the conda package: `vagrant ssh -c "source miniconda/bin/activate ~/miniconda; binstar -t $BINSTAR_TOKEN upload --force -u omnia -p openmm-dev $HOME/miniconda/conda-bld/linux-64/openmm-dev-*"`

# NOTES

* May need `vagrant plugin install vagrant-vbguest` after `vagrant up`.
