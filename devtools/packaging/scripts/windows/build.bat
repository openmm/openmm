REM "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64

cd openmm
mkdir build
cd build

set FFTW=C:\Program Files\Miniconda3\pkgs\fftw3f-3.3.4-vc14_2\Library
"C:\Program Files\CMake\bin\cmake.exe" .. -G "NMake Makefiles JOM" -DCMAKE_BUILD_TYPE=Release -DOPENMM_GENERATE_API_DOCS=ON -DFFTW_INCLUDES="%FFTW%/include" -DFFTW_LIBRARY="%FFTW%/lib/libfftw3f-3.lib"

jom
jom PythonInstall
jom C++ApiDocs
jom PythonApiDocs
# jom sphinxpdf
jom install
jom PythonBdist
