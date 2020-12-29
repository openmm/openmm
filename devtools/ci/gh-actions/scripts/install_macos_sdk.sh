# Install an older MacOS SDK
# This should guarantee OpenMM builds with extended compatibility across MacOS versions
# Adapted from conda-forge-ci-setup scripts:
# * https://github.com/conda-forge/conda-forge-ci-setup-feedstock/blob/dde296e/recipe/run_conda_forge_build_setup_osx
# * https://github.com/conda-forge/conda-forge-ci-setup-feedstock/blob/dde296e/recipe/download_osx_sdk.sh
#
# Some possible updates might involve upgrading the download link to future MacOS releases (10.15 to something else),
# depending on the version provided by the CI

OSX_SDK_DIR="$(xcode-select -p)/Platforms/MacOSX.platform/Developer/SDKs"
export MACOSX_DEPLOYMENT_TARGET=10.9
export MACOSX_SDK_VERSION=10.9

export CMAKE_OSX_SYSROOT="${OSX_SDK_DIR}/MacOSX${MACOSX_SDK_VERSION}.sdk"

if [[ ! -d ${CMAKE_OSX_SYSROOT}} ]]; then
    echo "Downloading ${MACOSX_SDK_VERSION} sdk"
    curl -L -O --connect-timeout 5 --max-time 10 --retry 5  --retry-delay 0 --retry-max-time 40 --retry-connrefused --retry-all-errors \
        https://github.com/phracker/MacOSX-SDKs/releases/download/10.15/MacOSX${MACOSX_SDK_VERSION}.sdk.tar.xz
    tar -xf MacOSX${MACOSX_SDK_VERSION}.sdk.tar.xz -C "$(dirname ${CMAKE_OSX_SYSROOT})"
fi

if [[ "$MACOSX_DEPLOYMENT_TARGET" == 10.* ]]; then
# set minimum sdk version to our target
plutil -replace MinimumSDKVersion -string ${MACOSX_SDK_VERSION} $(xcode-select -p)/Platforms/MacOSX.platform/Info.plist
plutil -replace DTSDKName -string macosx${MACOSX_SDK_VERSION}internal $(xcode-select -p)/Platforms/MacOSX.platform/Info.plist
fi

echo "MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}" >> ${GITHUB_ENV}
echo "CMAKE_OSX_SYSROOT=${MACOSX_DEPLOYMENT_TARGET}" >> ${GITHUB_ENV}
echo "CMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}" >> ${GITHUB_ENV}