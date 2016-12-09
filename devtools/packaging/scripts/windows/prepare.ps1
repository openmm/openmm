cd C:\Users\vagrant

# Install everything we can with choco.

choco install -y doxygen.portable swig cmake doxygen.install vcbuildtools git jom

# Install CUDA.

wget https://developer.nvidia.com/compute/cuda/8.0/prod/network_installers/cuda_8.0.44_win10_network-exe -UseBasicParsing -OutFile cuda_8.0.44_win10_network.exe
.\cuda_8.0.44_win10_network.exe -s compiler_8.0 cudart_8.0 cufft_8.0 cufft_dev_8.0 nvrtc_8.0 nvrtc_dev_8.0 | Out-Null

# Install Miniconda.

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe -UseBasicParsing -OutFile Miniconda3-latest-Windows-x86_64.exe
.\Miniconda3-latest-Windows-x86_64.exe /S /D=C:\Miniconda3 | Out-Null

# Install software with conda.

& "C:\Miniconda3\Scripts\conda.exe" install -y -c omnia fftw3f jinja2 lxml sphinx sphinxcontrib-autodoc_doxygen sphinxcontrib-lunrsearch conda-build
& "C:\Miniconda3\Scripts\pip.exe" install sphinxcontrib.bibtex