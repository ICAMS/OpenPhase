cmake_minimum_required(VERSION 3.10)

project(mpi_wrapper
    #VERSION 1.2.0.0
    DESCRIPTION "mpi_wrapper is custom interface for MPI and FFTW libraries."
    #    HOMEPAGE_URL "https://www.openphase.rub.de"
    LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(MPI REQUIRED)

if(WIN32)
    find_package(MKL REQUIRED)
    set(MKL_FFTW_INTERFACE  ${MKL_ROOT}/include/fftw)
    message(STATUS "Include MKL FFTW Interface: ${MKL_FFTW_INTERFACE}")
else()
    find_package(fftw REQUIRED)
endif()

# libLIB_WRAPPER
set(LIB_WRAPPER mpi_wrapper)
add_library(${LIB_WRAPPER} SHARED mpi_wrapper.cpp)

# libLIB_WRAPPER
if(NOT WIN32)
    target_include_directories(${LIB_WRAPPER} PUBLIC ${FFTW_INCLUDE_DIR})
    target_link_libraries(${LIB_WRAPPER} PUBLIC fftw3 fftw3_mpi)
endif()
target_link_libraries(${LIB_WRAPPER} PUBLIC MPI::MPI_CXX)

# Installing
#------------------------------------------------------------------------------
install(TARGETS ${LIB_WRAPPER} LIBRARY)

# Packing
#------------------------------------------------------------------------------
set(CPACK_PACKAGE_VENDOR "ICAMS")
include(CPack)
