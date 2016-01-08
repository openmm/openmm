:: Use python version to select which Visual Studio to use
:: For win-64, we'll need more, since those are separate compilers
:: Build in subdirectory.
mkdir build
cd build

set CMAKE_FLAGS=-DCMAKE_INSTALL_PREFIX=%PREFIX%
set CMAKE_FLAGS=%CMAKE_FLAGS% -DOPENMM_BUILD_PME_PLUGIN=ON
set CMAKE_FLAGS=%CMAKE_FLAGS% -DFFTW_LIBRARY=%LIBRARY_LIB%\libfftwf-3.3.lib
set CMAKE_FLAGS=%CMAKE_FLAGS% -DFFTW_INCLUDES=%LIBRARY_INC%
set CMAKE_FLAGS=%CMAKE_FLAGS% -DCMAKE_BUILD_TYPE=Release
set CMAKE_FLAGS=%CMAKE_FLAGS% -DOPENCL_INCLUDE_DIR="C:/Program Files (x86)/AMD APP SDK/2.9-1/include"
set CMAKE_FLAGS=%CMAKE_FLAGS% -DOPENCL_LIBRARY="C:/Program Files (x86)/AMD APP SDK/2.9-1/lib/x86_64/OpenCL.lib"

cmake -G "NMake Makefiles" %CMAKE_FLAGS% ..

jom all DoxygenApiDocs :: sphinxpdf
jom install
if errorlevel 1 exit 1


set OPENMM_INCLUDE_PATH=%PREFIX%\include
set OPENMM_LIB_PATH=%PREFIX%\lib
cd python
%PYTHON% setup.py install
cd ..

:: Put examples into an appropriate subdirectory.
mkdir %PREFIX%\share\openmm
move %PREFIX%\examples %PREFIX%\share\openmm

:: Put docs into a subdirectory.
cd %PREFIX%\docs
mkdir openmm
move *.html openmm 
:: move *.pdf openmm
:: move api-* openmm

if errorlevel 1 exit 1
