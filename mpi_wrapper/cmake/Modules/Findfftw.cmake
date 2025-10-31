# Find FFTW
# Try to find the FFTW libraries

find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARY fftw3)

if (NOT FFTW_INCLUDE_DIR)
  message(FATAL_ERROR "FFTW include directory not found.")
else()
  message(STATUS "Found FFTW include directory: ${FFTW_INCLUDE_DIR}")
endif()

if (NOT FFTW_LIBRARY)
  message(FATAL_ERROR "FFTW library not found.")
else()
  message(STATUS "Found FFTW library: ${FFTW_LIBRARY}")
endif()


# If found, set the FFTW include and library variables
set(FFTW_FOUND TRUE)
set(FFTW_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
set(FFTW_LIBRARY ${FFTW_LIBRARY})
