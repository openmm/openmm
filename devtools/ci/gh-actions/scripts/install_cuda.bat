:: This script installs CUDA on Windows.
:: It downloads the offline installer from the Nvidia servers
:: and the relevant patches, if applicable.
:: It uses the default installation path, which is exported as CUDA_PATH
:: For CMake compatibility, CUDA_TOOLKIT_ROOT_DIR is also exported
:: It expects a %CUDA_VERSION% environment variable, set to major.minor (e.g. 10.0)

:: We define a default subset of components to be installed for faster installation times
:: and reduced storage usage (CI is limited to 10GB). Full list of components is available at
:: https://docs.nvidia.com/cuda/archive/%CUDA_VERSION%/cuda-installation-guide-microsoft-windows/index.html
set "VAR=nvcc_%CUDA_VERSION% cuobjdump_%CUDA_VERSION% nvprune_%CUDA_VERSION% cupti_%CUDA_VERSION%"
set "VAR=%VAR% memcheck_%CUDA_VERSION% nvdisasm_%CUDA_VERSION% nvprof_%CUDA_VERSION% cublas_%CUDA_VERSION%"
set "VAR=%VAR% cublas_dev_%CUDA_VERSION% cudart_%CUDA_VERSION% cufft_%CUDA_VERSION% cufft_dev_%CUDA_VERSION%"
set "VAR=%VAR% curand_%CUDA_VERSION% curand_dev_%CUDA_VERSION% cusolver_%CUDA_VERSION% cusolver_dev_%CUDA_VERSION%"
set "VAR=%VAR% cusparse_%CUDA_VERSION% cusparse_dev_%CUDA_VERSION% npp_%CUDA_VERSION% npp_dev_%CUDA_VERSION%"
set "VAR=%VAR% nvrtc_%CUDA_VERSION% nvrtc_dev_%CUDA_VERSION% nvml_dev_%CUDA_VERSION%"
set "VAR=%VAR% visual_studio_integration_%CUDA_VERSION%"
set "CUDA_COMPONENTS=%VAR%"

if "%CUDA_VERSION%" == "9.2"  goto cuda92
if "%CUDA_VERSION%" == "10.0" goto cuda100
if "%CUDA_VERSION%" == "10.1" goto cuda101
if "%CUDA_VERSION%" == "10.2" goto cuda102
if "%CUDA_VERSION%" == "11.0" goto cuda110
if "%CUDA_VERSION%" == "11.1" goto cuda111
if "%CUDA_VERSION%" == "11.2" goto cuda112

echo CUDA '%CUDA_VERSION%' is not supported
exit /b 1

:: Define URLs per version
:cuda92
set "CUDA_NETWORK_INSTALLER_URL=https://developer.nvidia.com/compute/cuda/9.2/Prod2/network_installers2/cuda_9.2.148_win10_network"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=2bf9ae67016867b68f361bf50d2b9e7b"
set "CUDA_INSTALLER_URL=https://developer.nvidia.com/compute/cuda/9.2/Prod2/local_installers2/cuda_9.2.148_win10"
set "CUDA_INSTALLER_CHECKSUM=f6c170a7452098461070dbba3e6e58f1"
set "CUDA_PATCH_URL=https://developer.nvidia.com/compute/cuda/9.2/Prod2/patches/1/cuda_9.2.148.1_windows"
set "CUDA_PATCH_CHECKSUM=09e20653f1346d2461a9f8f1a7178ba2"
set "CUDA_COMPONENTS=%CUDA_COMPONENTS% nvgraph_%CUDA_VERSION% nvgraph_dev_%CUDA_VERSION%"
goto cuda_common


:cuda100
set "CUDA_NETWORK_INSTALLER_URL=https://developer.nvidia.com/compute/cuda/10.0/Prod/network_installers/cuda_10.0.130_win10_network"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=3312deac9c939bd78d0e7555606c22fc"
set "CUDA_INSTALLER_URL=https://developer.nvidia.com/compute/cuda/10.0/Prod/local_installers/cuda_10.0.130_411.31_win10"
set "CUDA_INSTALLER_CHECKSUM=90fafdfe2167ac25432db95391ca954e"
set "CUDA_COMPONENTS=%CUDA_COMPONENTS% nvgraph_%CUDA_VERSION% nvgraph_dev_%CUDA_VERSION%"
goto cuda_common


:cuda101
set "CUDA_NETWORK_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/10.1/Prod/network_installers/cuda_10.1.243_win10_network.exe"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=fae0c958440511576691b825d4599e93"
set "CUDA_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/10.1/Prod/local_installers/cuda_10.1.243_426.00_win10.exe"
set "CUDA_INSTALLER_CHECKSUM=b54cf32683f93e787321dcc2e692ff69"
set "CUDA_COMPONENTS=%CUDA_COMPONENTS% nvgraph_%CUDA_VERSION% nvgraph_dev_%CUDA_VERSION%"
goto cuda_common


:cuda102
set "CUDA_NETWORK_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/10.2/Prod/network_installers/cuda_10.2.89_win10_network.exe"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=60e0f16845d731b690179606f385041e"
set "CUDA_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda_10.2.89_441.22_win10.exe"
set "CUDA_INSTALLER_CHECKSUM=d9f5b9f24c3d3fc456a3c789f9b43419"
set "CUDA_PATCH_URL=http://developer.download.nvidia.com/compute/cuda/10.2/Prod/patches/1/cuda_10.2.1_win10.exe"
set "CUDA_PATCH_CHECKSUM=9d751ae129963deb7202f1d85149c69d"
set "CUDA_COMPONENTS=%CUDA_COMPONENTS% nvgraph_%CUDA_VERSION% nvgraph_dev_%CUDA_VERSION%"
goto cuda_common


:cuda110
set "CUDA_NETWORK_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/11.0.3/network_installers/cuda_11.0.3_win10_network.exe"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=1b88bf7bb8e50207bbb53ed2033f93f3"
set "CUDA_INSTALLER_URL=http://developer.download.nvidia.com/compute/cuda/11.0.3/local_installers/cuda_11.0.3_451.82_win10.exe"
set "CUDA_INSTALLER_CHECKSUM=80ae0fdbe04759123f3cab81f2aadabd"
goto cuda_common


:cuda111
set "CUDA_NETWORK_INSTALLER_URL=https://developer.download.nvidia.com/compute/cuda/11.1.1/network_installers/cuda_11.1.1_win10_network.exe"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=7e36e50ee486a84612adfd85500a9971"
set "CUDA_INSTALLER_URL=https://developer.download.nvidia.com/compute/cuda/11.1.1/local_installers/cuda_11.1.1_456.81_win10.exe"
set "CUDA_INSTALLER_CHECKSUM=a89dfad35fc1adf02a848a9c06cfff15"
goto cuda_common


:cuda112
set "CUDA_NETWORK_INSTALLER_URL=https://developer.download.nvidia.com/compute/cuda/11.2.0/network_installers/cuda_11.2.0_win10_network.exe"
set "CUDA_NETWORK_INSTALLER_CHECKSUM=ab02a25eed1201cc3e414be943a242df"
set "CUDA_INSTALLER_URL=https://developer.download.nvidia.com/compute/cuda/11.2.0/local_installers/cuda_11.2.0_460.89_win10.exe"
set "CUDA_INSTALLER_CHECKSUM=92f38c37ce9c6c11d27c10701b040256"
goto cuda_common


:: The actual installation logic
:cuda_common

::We expect this CUDA_PATH
set "CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v%CUDA_VERSION%"

echo Downloading CUDA version %CUDA_VERSION% installer from %CUDA_INSTALLER_URL%
echo Expected MD5: %CUDA_INSTALLER_CHECKSUM%

:: Download installer
curl --retry 3 -k -L %CUDA_INSTALLER_URL% --output cuda_installer.exe
if errorlevel 1 (
    echo Problem downloading installer...
    exit /b 1
)
:: Check md5
openssl md5 cuda_installer.exe | findstr %CUDA_INSTALLER_CHECKSUM%
if errorlevel 1 (
    echo Checksum does not match!
    exit /b 1
)
:: Run installer
start /wait cuda_installer.exe -s %CUDA_COMPONENTS%
if errorlevel 1 (
    echo Problem installing CUDA toolkit...
    exit /b 1
)
del cuda_installer.exe

:: If patches are needed, download and apply
if not "%CUDA_PATCH_URL%"=="" (
    echo This version requires an additional patch
    curl --retry 3 -k -L %CUDA_PATCH_URL% --output cuda_patch.exe
    if errorlevel 1 (
        echo Problem downloading patch installer...
        exit /b 1
    )
    openssl md5 cuda_patch.exe | findstr %CUDA_PATCH_CHECKSUM%
    if errorlevel 1 (
        echo Checksum does not match!
        exit /b 1
    )
    start /wait cuda_patch.exe -s
    if errorlevel 1 (
        echo Problem running patch installer...
        exit /b 1
    )
    del cuda_patch.exe
)

:: This should exist by now!
if not exist "%CUDA_PATH%\bin\nvcc.exe" (
    echo CUDA toolkit installation failed!
    exit /b 1
)

echo CUDA_PATH=%CUDA_PATH% >> %GITHUB_ENV%
echo CUDA_TOOLKIT_ROOT_DIR=%CUDA_PATH:\=/% >> %GITHUB_ENV%

:: Notes about nvcuda.dll
:: ----------------------
:: We should also provide the drivers (nvcuda.dll), but the installer will not
:: proceed without a physical Nvidia card attached (not the case in the CI).
:: Expanding `<installer.exe>\Display.Driver\nvcuda.64.dl_` to `C:\Windows\System32`
:: does not work anymore (.dl_ files are not PE-COFF according to Dependencies.exe).
:: Forcing this results in a DLL error 193. Basically, there's no way to provide
:: ncvuda.dll in a GPU-less machine without breaking the EULA (aka zipping nvcuda.dll
:: from a working installation).