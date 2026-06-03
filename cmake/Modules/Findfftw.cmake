# Find FFTW
# Try to find the FFTW libraries

set(FFTW_HINT_INCLUDE_DIRS)
set(FFTW_HINT_LIBRARY_DIRS)

if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND FFTW_HINT_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/include")
  list(APPEND FFTW_HINT_LIBRARY_DIRS "$ENV{CONDA_PREFIX}/lib")
endif()

find_path(FFTW_INCLUDE_DIR fftw3.h
  HINTS ${FFTW_HINT_INCLUDE_DIRS}
)
find_library(FFTW_LIBRARY fftw3
  HINTS ${FFTW_HINT_LIBRARY_DIRS}
)

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
