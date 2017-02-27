# Packaging OpenMM into ZIP installers

## Source

Start the docker container:
```bash
docker run -i -t --rm -v `pwd`:/io jchodera/omnia-build-box:cuda80-amd30-clang38 bash
```
Inside the docker container:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/source/prepare.sh       
source openmm/devtools/packaging/scripts/source/build.sh
source openmm/devtools/packaging/scripts/source/package.sh
# Recover the packages to host directory
cp packaging/compressed/* /io
```

## Linux

Start the docker container:
```bash
docker run -i -t --rm -v `pwd`:/io jchodera/omnia-build-box:cuda80-amd30-clang38 bash
```
Inside the docker container:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/linux/prepare.sh
source openmm/devtools/packaging/scripts/linux/build.sh
source openmm/devtools/packaging/scripts/linux/package.sh
# Recover the packages to host directory
cp packaging/compressed/* /io
```

## OS X

On an `osx` machine with XCode and the OS X 10.9 frameworks installed:
```bash
# Clone the OpenMM beta or release candidate tag $TAG
git clone https://github.com/pandegroup/openmm.git
cd openmm; git checkout $TAG; cd ..
# Build and package
source openmm/devtools/packaging/scripts/osx/prepare.sh
source openmm/devtools/packaging/scripts/osx/build.sh
source openmm/devtools/packaging/scripts/osx/package.sh
```

