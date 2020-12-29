# This script installs AMD's SDK 3.0 to provide their OpenCL implementation
# * Installation path will be ${GITHUB_WORKSPACE}/AMDAPPSDK

set -euxo pipefail


wget -q --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --tries 5 \
    http://s3.amazonaws.com/omnia-ci/AMD-APP-SDKInstaller-v3.0.130.135-GA-linux64.tar.bz2
tar -xjf AMD-APP-SDK*.tar.bz2

AMDAPPSDK=${GITHUB_WORKSPACE}/AMDAPPSDK
export OPENCL_VENDOR_PATH=${AMDAPPSDK}/etc/OpenCL/vendors

mkdir -p ${OPENCL_VENDOR_PATH}
sh AMD-APP-SDK*.sh --tar -xf -C ${AMDAPPSDK}
echo libamdocl64.so > ${OPENCL_VENDOR_PATH}/amdocl64.icd

export LD_LIBRARY_PATH=${AMDAPPSDK}/lib/x86_64:${LD_LIBRARY_PATH:-}
chmod +x ${AMDAPPSDK}/bin/x86_64/clinfo
${AMDAPPSDK}/bin/x86_64/clinfo
sudo apt-get install -y libgl1-mesa-dev

echo "OPENCL_VENDOR_PATH=${OPENCL_VENDOR_PATH}" >> ${GITHUB_ENV}
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> ${GITHUB_ENV}