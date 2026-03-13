set(CMAKE_C_COMPILER icx)
set(CMAKE_CXX_COMPILER icx)

# optional but VERY recommended
set(CMAKE_LINKER lld-link)

# Windows runtime compatibility
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
