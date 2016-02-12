#!/bin/bash

# Packaging script for Linux distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

# CONFIGURE HERE
export PACKAGE_DIR="packaging" # directory to stuff packaged source distribution
export VERSION=$(sed -nr "s/OPENMM_VERSION:STRING=(.*)/\1/p" build/CMakeCache.txt)
export PACKAGE_SUBDIR="OpenMM-${VERSION}-Linux" # directory where distribution will be unpacked
export DISTRO_PREFIX="OpenMM-${VERSION}-Linux" # prefix for source distribution (e.g. ${DISTRIBUTION_NAME}.zip)

# Perform all work in a work directory.
cd work

# Clean up.
rm -rf $PACKAGE_DIR

# Make a directory to contain packaged source distribution
mkdir $PACKAGE_DIR
mkdir $PACKAGE_DIR/$PACKAGE_SUBDIR
for filename in $( cat openmm/devtools/packaging/manifests/binary/manifest.txt ); do
   CMD="cp -r install/$filename $PACKAGE_DIR/$PACKAGE_SUBDIR"
   echo $CMD
   `$CMD`
done

# Add the install.sh script
CMD="cp -r openmm/devtools/packaging/install.sh $PACKAGE_DIR/$PACKAGE_SUBDIR"
echo $CMD
`$CMD`

# Make Python source distribution.
echo "Building Python source distribution..."
pushd .
cd build
make PythonSdist
cd python/dist
tar zxf OpenMM-${VERSION}.tar.gz
mv OpenMM-${VERSION} python
popd
cp -r build/python/dist/python $PACKAGE_DIR/$PACKAGE_SUBDIR

# Create archives.
cd $PACKAGE_DIR
mkdir compressed
tar zcf compressed/${DISTRO_PREFIX}.tgz $PACKAGE_SUBDIR
zip -r compressed/${DISTRO_PREFIX}.zip $PACKAGE_SUBDIR
cd ..
