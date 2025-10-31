set(app_name Flame0D)
find_package(Threads)
find_package(MPI REQUIRED)

include_directories("$ENV{CANTERA_PATH}/include")
link_directories("$ENV{CANTERA_PATH}/lib")

add_executable(${app_name} ${app_name}.cpp)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ProjectInput.opi DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/SettingsInput.opi DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# If the example calls MPI_Init/Finalize directly, link MPI too:
# target_link_libraries(${app_name} PRIVATE MPI::MPI_CXX)

# Link OpenPhase (this already propagates its PUBLIC include dirs/libs)
target_link_libraries(${app_name} PRIVATE OpenPhase)

# Only add direct deps the example itself uses explicitly
target_link_libraries(${app_name} PRIVATE fftw3 cantera Threads::Threads)

# Ensure C++17 on the example (matches the library)
target_compile_features(${app_name} PRIVATE cxx_std_17)
