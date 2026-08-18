cmake_minimum_required(VERSION 3.15)
project(gas_transmittance LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)

if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

pybind11_add_module(gas_transmittance
    src/gas_transmittance.cpp
    src/gas_transmittance_impl.cpp
)

target_include_directories(gas_transmittance PRIVATE include)
target_link_libraries(gas_transmittance PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS gas_transmittance DESTINATION gas_corrections_lib/gas_corrections_cpp)