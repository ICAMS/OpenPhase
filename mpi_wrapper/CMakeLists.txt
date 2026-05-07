cmake_minimum_required(VERSION 3.10)

project(mpi_wrapper
    DESCRIPTION "MPI + FFTW wrapper"
    LANGUAGES CXX C
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(MPI REQUIRED)

add_library(mpi_wrapper SHARED mpi_wrapper.cpp)

if(WIN32)
    # Integer interface (HPC default)
    set(MKL_INTERFACE "ilp64" CACHE STRING "MKL integer interface (lp64 or ilp64)")

    # Set MKL integer interface (default: ilp64 for HPC)
    if(MKL_INTERFACE STREQUAL "ilp64")
        target_compile_definitions(mpi_wrapper PUBLIC MKL_ILP64)
    elseif(MKL_INTERFACE STREQUAL "lp64")
        # LP64 explicitly requested
    else()
        message(FATAL_ERROR "MKL_INTERFACE must be 'lp64' or 'ilp64'")
    endif()

    cmake_policy(SET CMP0074 NEW)

    find_package(MKL REQUIRED)

    # FFTW headers from MKL
    target_include_directories(mpi_wrapper PUBLIC
        ${MKL_ROOT}/include
        ${MKL_ROOT}/include/fftw
    )

    if(MKL_INTERFACE STREQUAL "ilp64")
        set(FFTW_MPI_LIB ${MKL_ROOT}/lib/mkl_blacs_intelmpi_ilp64.lib)
    else()
        set(FFTW_MPI_LIB ${MKL_ROOT}/lib/mkl_blacs_intelmpi_lp64.lib)
    endif()

    target_link_libraries(mpi_wrapper PUBLIC
        MPI::MPI_CXX
        ${FFTW_MPI_LIB}
        ${MKL_LIBRARIES}
    )
else()
    find_package(fftw REQUIRED)
    target_include_directories(mpi_wrapper PUBLIC ${FFTW_INCLUDE_DIR})
    target_link_libraries(mpi_wrapper PUBLIC MPI::MPI_CXX fftw3 fftw3_mpi)
endif()

install(TARGETS mpi_wrapper LIBRARY)
