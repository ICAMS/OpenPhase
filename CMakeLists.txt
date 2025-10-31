cmake_minimum_required(VERSION 3.10)

project(OpenPhase
    VERSION 2.0.0.0
    DESCRIPTION "OpenPhase is an open source software library for the phase-field simulations of complex scientific problems."
    #    HOMEPAGE_URL "https://www.openphase.rub.de"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(ENABLE_MPI "Enable MPI support" ON)
option(ENABLE_SENTINEL "Enable Sentinel (copy protection)" OFF)
option(ENABLE_ACADEMIC "Enable compilation of academic library version" ON)
option(ENABLE_BENCHMARKS "Enable compilation of benchmarks" ON)
option(ENABLE_EXAMPLES "Enable compilation of examples" ON)
option(ENABLE_DYNAMIC_LIBS "Enable compilation shared OpenPhase library" ON)
option(ENABLE_DYNAMIC_LINKING "Enable shared linking of dependencies" ON)

if (NOT ENABLE_DYNAMIC_LINKING AND ENABLE_OPENMP AND ENABLE_SENTINEL)
    message(FATAL_ERROR "Static linking of OpenMP and Sentinel is not supported.")
endif()

# Generate BuildInfo.h
set(BUILD_INFO_HEADER "${PROJECT_SOURCE_DIR}/include/BuildInfo.h")

# Get the current Git commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE GIT_REV_PARSE_RESULT
    OUTPUT_VARIABLE GIT_COMMIT_SHA
    ERROR_QUIET
)
string(STRIP "${GIT_COMMIT_SHA}" GIT_COMMIT_SHA)  # Remove any trailing newlines or spaces

# Get the current build time
execute_process(
    COMMAND date
    OUTPUT_VARIABLE BUILD_TIME
    ERROR_QUIET
)
string(STRIP "${BUILD_TIME}" BUILD_TIME)  # Remove any trailing newlines or spaces

# Create the BuildInfo.h file
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/BuildInfo.h.in"
    "${BUILD_INFO_HEADER}"
    @ONLY
)

# Find dependencies
#------------------------------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

# FFTW3
if(WIN32)
    find_package(MKL REQUIRED)
    set(MKL_FFTW_INTERFACE  ${MKL_ROOT}/include/fftw)
    message(STATUS "MKL_FFTW_INTERFACE: ${MKL_FFTW_INTERFACE}")
    set(FFTW3_LIBRARIES ${MKL_LIBRARIES})
    message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qmkl")
else()
    if (ENABLE_DYNAMIC_LINKING)
        find_package(fftw REQUIRED)
        set(FFTW3_LIBRARIES fftw3)
    else()
        find_library(FFTW3_LIBRARIES libfftw3.a)
        message(STATUS "FOUND fftw3: ${FFTW3_LIBRARIES}")
    endif()
endif()

# OpenMP
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        if(ENABLE_DYNAMIC_LINKING)
            set(OPENMP_LIBRARIES OpenMP::OpenMP_CXX)
            if(NOT WIN32)
                list(APPEND FFTW3_LIBRARIES fftw3_omp)
            endif()
        else()
            if(NOT WIN32)
                find_library(OPENMP_LIBRARIES libgomp.a HINTS /usr/lib/gcc/x86_64-linux-gnu/*/)
                find_library(FFTW3_OMP_LIBRARIES libfftw3_omp.a)
                list(APPEND FFTW3_LIBRARIES ${FFTW3_OMP_LIBRARIES})
            endif()
            message(STATUS "FOUND OpenMP: ${OPENMP_LIBRARIES}")
        endif()
    endif()
endif()
if (NOT OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_OPENMP")
endif()

# MPI
if (ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        if(ENABLE_DYNAMIC_LINKING)
            set(MPI_LIBRARIES MPI::MPI_CXX)
            if(NOT Win32)
                list(APPEND FFTW3_LIBRARIES fftw3_mpi)
            endif()
        else()
            message(FATAL_ERROR "Static linking of MPI is not yet supported.")
        endif()
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DMPI_PARALLEL")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMPI_PARALLEL")
    endif()
endif()

# Sentinel
if (ENABLE_SENTINEL)
    if(WIN32)
        set(SENTINEL_LIBRARY ${PROJECT_SOURCE_DIR}/external/Sentinel/C/x86_64/hasp_windows_x64_96982.lib)
    else()
        set(SENTINEL_LIBRARY ${PROJECT_SOURCE_DIR}/external/Sentinel/C/x86_64/libhasp_linux_x86_64_96982.a)
    endif()
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOSENTINEL")
endif()


# Cantera
if(DEFINED ENV{CANTERA_PATH})
    set(CANTERA_FOUND TRUE)
    message(STATUS "FOUND Cantera")
    include_directories("$ENV{CANTERA_PATH}/include")
    link_directories("$ENV{CANTERA_PATH}/lib")
    find_package(Threads REQUIRED)
    find_package(fmt REQUIRED)
    set(CANTERA_LIBRARIES  cantera Threads::Threads)
    message(STATUS "CANTERA_LIBRARIES: ${CANTERA_LIBRARIES}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DCANTERA")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DCANTERA")
else()
    set(CANTERA_FOUND FALSE)
    #message(STATUS "Cantera path not found.")
endif()

# Set flags
#------------------------------------------------------------------------------

# define build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug) # Set debug mode as standard for development branch!
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -DDEBUG")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (ENABLE_ACADEMIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DACADEMIC")
endif()


# Includes
#------------------------------------------------------------------------------
set(OP_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external
    ${PROJECT_SOURCE_DIR}/mpi_wrapper
    ${FFTW_INCLUDE_DIR})

if(WIN32)
    list(APPEND OP_INCLUDES ${MKL_FFTW_INTERFACE})
endif()

if (NOT ENABLE_ACADEMIC)
    list(APPEND OP_INCLUDES ${PROJECT_SOURCE_DIR}/include_closed)
endif()

# Sources
#------------------------------------------------------------------------------
file(GLOB OP_SOURCES
    ${PROJECT_SOURCE_DIR}/external/*.c
    ${PROJECT_SOURCE_DIR}/external/*.cpp
    ${PROJECT_SOURCE_DIR}/external/*/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*/*/*.cpp)
list(FILTER OP_SOURCES EXCLUDE REGEX "DislocationModel.cpp$")
list(FILTER OP_SOURCES EXCLUDE REGEX "DislocationModelLine.cpp$")
list(FILTER OP_SOURCES EXCLUDE REGEX "DoubleObstacle_SPF.cpp$")
list(FILTER OP_SOURCES EXCLUDE REGEX "FiniteInterfaceDissipation.cpp$")
list(FILTER OP_SOURCES EXCLUDE REGEX "Stoichiometry.cpp$")

if (NOT ENABLE_ACADEMIC)
    file(GLOB OP_SOURCES_CLOSED
        ${PROJECT_SOURCE_DIR}/src_closed/*.cpp
        ${PROJECT_SOURCE_DIR}/src_closed/*/*.cpp
        ${PROJECT_SOURCE_DIR}/src_closed/*/*/*.cpp)
    if(WIN32)# OR INCLUDE_CLOSED_SOURCES)
        message(STATUS "Compile closed sources into ${PROJECT_NAME} library")
        list(APPEND OP_SOURCES ${OP_SOURCES_CLOSED})
    endif()
endif()

#message(STATUS "OP_SOURCES: ${OP_SOURCES}")
#message(STATUS "OP_SOURCES_CLOSED: ${OP_SOURCES_CLOSED}")

# Targets
#------------------------------------------------------------------------------

set(OP_LIBS)

# libLIB_WRAPPER
if (ENABLE_MPI AND MPI_FOUND)
    set(LIB_WRAPPER mpi_wrapper)
    add_subdirectory(mpi_wrapper)
    list(APPEND OP_LIBS ${LIB_WRAPPER})
endif()

# libOpenPhase
set(LIB_OPENPHASE OpenPhase)
if (ENABLE_DYNAMIC_LIBS)
    add_library(${LIB_OPENPHASE} SHARED ${OP_SOURCES})
else()
    add_library(${LIB_OPENPHASE} STATIC ${OP_SOURCES})
endif()
target_include_directories(${LIB_OPENPHASE} PUBLIC ${OP_INCLUDES})
set_target_properties(${LIB_OPENPHASE} PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_definitions(OpenPhase PRIVATE OPENPHASE_EXPORTS)
list(APPEND OP_LIBS ${LIB_OPENPHASE})


# libOpenPhaseClosed
if (NOT WIN32 AND NOT ENABLE_ACADEMIC)
    set(LIB_OPENPHASE_CLOSED OpenPhaseClosed)
    add_library(${LIB_OPENPHASE_CLOSED} SHARED ${OP_SOURCES_CLOSED})
    target_include_directories(${LIB_OPENPHASE_CLOSED} PUBLIC ${OP_INCLUDES})
    set_target_properties(${LIB_OPENPHASE_CLOSED} PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
    list(APPEND OP_LIBS ${LIB_OPENPHASE_CLOSED})
endif()

# Link Libraries
#------------------------------------------------------------------------------

# libOpenPhase
target_link_libraries(${LIB_OPENPHASE} PUBLIC ${FFTW3_LIBRARIES})
target_link_libraries(${LIB_OPENPHASE} PUBLIC ${OPENMP_LIBRARIES})
target_link_libraries(${LIB_OPENPHASE} PUBLIC ${LIB_WRAPPER})
if (CANTERA_FOUND)
    target_link_libraries(${LIB_OPENPHASE} PUBLIC ${CANTERA_LIBRARIES})
endif()

if (NOT WIN32 AND NOT ENABLE_ACADEMIC)
    target_link_libraries(${LIB_OPENPHASE} PUBLIC ${LIB_OPENPHASE_CLOSED})
    target_link_libraries(${LIB_OPENPHASE_CLOSED} PUBLIC ${FFTW3_LIBRARIES})
    target_link_libraries(${LIB_OPENPHASE_CLOSED} PUBLIC ${OPENMP_LIBRARIES})
    target_link_libraries(${LIB_OPENPHASE_CLOSED} PUBLIC ${SENTINEL_LIBRARY})
endif()

# Testing
#------------------------------------------------------------------------------
enable_testing()

# Subdirectories (Executables)
#------------------------------------------------------------------------------

# define common operation for all executables of OpenPhase
function(add_openphase_executable EXENAME)
    add_executable(${EXENAME} ${ARGN})
    target_link_libraries(${EXENAME} PRIVATE ${OP_LIBS})

    # ensure OpenPhase DLL is built first
    add_dependencies(${EXENAME} ${LIB_OPENPHASE})

    # Copy opi files
    file(GLOB OPI_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.opi")
    if (OPI_FILES)
        foreach(OPI_FILE IN LISTS OPI_FILES)
            get_filename_component(OPI_NAME ${OPI_FILE} NAME)
            configure_file(
                ${OPI_FILE}
                ${CMAKE_CURRENT_BINARY_DIR}/${OPI_NAME}
                COPYONLY
            )
        endforeach()
    endif()

    # Copy DLL on Windows
    if (WIN32)
        add_custom_command(TARGET ${EXENAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${LIB_OPENPHASE}>"
                $<TARGET_FILE_DIR:${EXENAME}>
        )
    endif()
endfunction()

if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

# Installing
#------------------------------------------------------------------------------
install(TARGETS ${LIB_OPENPHASE} LIBRARY)
#install(TARGETS ${LIB_OPENPHASE_CLOSED} LIBRARY)

# Packing
#------------------------------------------------------------------------------
set(CPACK_PACKAGE_VENDOR "ICAMS")
include(CPack)
