set(app_name PoiseFlow)

find_package(MPI REQUIRED)

add_executable(${app_name} ${app_name}.cpp)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ProjectInput.opi DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(${app_name} PUBLIC OpenPhase)
