# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES        - where to find fftw3.h
#  FFTW_LIBRARY         - the main FFTW library.
#  FFTW_THREADS_LIBRARY - the FFTW multithreading support library.
#  FFTW_FOUND           - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

find_library (FFTW_LIBRARY NAMES fftw3f)
find_library (FFTW_THREADS_LIBRARY NAMES fftw3f_threads)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARY FFTW_THREADS_LIBRARY FFTW_INCLUDES)
