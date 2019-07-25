mkdir build
cd build

set FFTW=C:\Miniconda3\pkgs\fftw3f-3.3.4-vc14_2\Library
set APPSDK=C:\Program Files (x86)\AMD APP SDK\2.9-1
"C:\Program Files\CMake\bin\cmake.exe" .. -G "NMake Makefiles JOM" -DCMAKE_BUILD_TYPE=Release -DOPENMM_GENERATE_API_DOCS=ON ^
    -DOPENCL_INCLUDE_DIR="%APPSDK%\include" -DOPENCL_LIBRARY="%APPSDK%\lib\x86_64\OpenCL.lib" ^
    -DFFTW_INCLUDES="%FFTW%/include" -DFFTW_LIBRARY="%FFTW%/lib/libfftw3f-3.lib"

jom
jom PythonInstall
jom C++ApiDocs
jom PythonApiDocs
REM jom sphinxpdf
jom install
jom PythonBdist
