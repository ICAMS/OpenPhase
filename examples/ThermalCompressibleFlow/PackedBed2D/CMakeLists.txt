set(app_name PB2D)

find_package(MPI REQUIRED)

add_executable(${app_name} ${app_name}.cpp)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ProjectInput.opi DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cylindersData.txt    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(${app_name} PUBLIC OpenPhase)
