set(app_name InterfaceDiffusionVolumeConservation)
add_openphase_executable(${app_name} ${app_name}.cpp)
add_test(${app_name} ${app_name})
