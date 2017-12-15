cd C:\Users\vagrant

# Install CUDA.

wget https://developer.nvidia.com/compute/cuda/9.0/Prod/network_installers/cuda_9.0.176_win10_network-exe -UseBasicParsing -OutFile cuda_9.0.176_win10_network.exe
.\cuda_9.0.176_win10_network.exe -s compiler_9.0 cudart_9.0 cufft_9.0 cufft_dev_9.0 nvrtc_9.0 nvrtc_dev_9.0 | Out-Null

# Install AMD APP SDK.

wget https://s3.amazonaws.com/omnia-ci/AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe -UseBasicParsing -OutFile AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe
.\AMD-APP-SDK-v2.9-1.599.381-GA-Full-windows-64.exe /S /v/qn | Out-Null

# Install Miniconda.

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe -UseBasicParsing -OutFile Miniconda3-latest-Windows-x86_64.exe
.\Miniconda3-latest-Windows-x86_64.exe /S /D=C:\Miniconda3 | Out-Null

# Install software with conda.

& "C:\Miniconda3\Scripts\conda.exe" install -y -c omnia fftw3f jinja2 lxml sphinx sphinxcontrib-autodoc_doxygen sphinxcontrib-lunrsearch conda-build anaconda-client
& "C:\Miniconda3\Scripts\pip.exe" install sphinxcontrib.bibtex

# Install software with choco.

choco install -y doxygen.portable swig cmake doxygen.install vcbuildtools git jom
