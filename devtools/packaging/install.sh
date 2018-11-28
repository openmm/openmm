#!/bin/sh

cd $(dirname $0)

# Ask the user for the install location and Python executable.

defaultInstallDir=/usr/local/openmm
printf "Enter install location (default=${defaultInstallDir}): "
read installDir
if [ -z ${installDir} ]
then
  installDir=${defaultInstallDir}
fi
defaultPythonBin=$(which python)
printf "Enter path to Python executable"
if [ ${defaultPythonBin} ]
then
  printf " (default=${defaultPythonBin})"
fi
printf ": "
read pythonBin
if [ -z ${pythonBin} ]
then
  pythonBin=${defaultPythonBin}
fi

# Make sure it's a supported Python version.

pythonOk=$(${pythonBin} -c "import sys; v=sys.version_info; print((v[0]==2 and v[1]>6) or v[0]>2)")
if [ ${pythonOk} != "True" ]
then
  echo "Unsupported Python version.  Only versions 2.7 and higher are supported."
  exit
fi

# Copy the files into place.

cp -R docs ${installDir}
cp -R include ${installDir}
cp -R lib ${installDir}
cp -R licenses ${installDir}

# Run the Python installer.

cd python
export OPENMM_INCLUDE_PATH=${installDir}/include
export OPENMM_LIB_PATH=${installDir}/lib
printenv
if ${pythonBin} setup.py build && ${pythonBin} setup.py install $@
then
    # Print instructions to the user.

    echo
    echo "Installation is complete.  You should now test your installation to make sure"
    echo "it is working correctly by typing the following command:"
    echo
    echo "python -m simtk.testInstallation"
else
    echo
    echo "INSTALLATION FAILED"
    echo
    echo "An error prevented the installation from completing.  See above for details." 
fi
