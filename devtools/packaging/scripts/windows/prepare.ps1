cd C:\Users\vagrant

# Install CUDA.

wget https://developer.nvidia.com/compute/cuda/10.1/Prod/network_installers/cuda_10.1.168_win10_network.exe -UseBasicParsing -OutFile cuda_10.1.168_win10_network.exe
.\cuda_10.1.168_win10_network.exe -s nvcc_10.1 cudart_10.1 cufft_10.1 cufft_dev_10.1 nvrtc_10.1 nvrtc_dev_10.1 | Out-Null

# Install AMD APP SDK.

wget https://s3.amazonaws.com/omnia-ci/AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe -UseBasicParsing -OutFile AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe
.\AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe /S /v/qn | Out-Null

# Install Miniconda.

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe -UseBasicParsing -OutFile Miniconda3-latest-Windows-x86_64.exe
.\Miniconda3-latest-Windows-x86_64.exe /S /D=C:\Miniconda3 | Out-Null
[Environment]::SetEnvironmentVariable("Path", $env:Path + ";C:\Miniconda3;C:\Miniconda3\Scripts;C:\Miniconda3\Library\bin", [EnvironmentVariableTarget]::User)

# Install software with conda.

& "C:\Miniconda3\Scripts\conda.exe" config --add channels omnia --add channels conda-forge
& "C:\Miniconda3\Scripts\conda.exe" install -y fftw3f==3.3.4=vc14_2 jinja2 lxml sphinx sphinxcontrib-autodoc_doxygen sphinxcontrib-lunrsearch conda-build anaconda-client
& "C:\Miniconda3\Scripts\pip.exe" install sphinxcontrib.bibtex

# Install software with choco.

choco install -y doxygen.portable swig cmake doxygen.install vcbuildtools git jom patch
