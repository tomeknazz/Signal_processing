cmake_minimum_required(VERSION 3.4...3.18)
project(Signal)

add_subdirectory(pybind11)
pybind11_add_module(Signal Signal_processing.cpp)
add_subdirectory(matplotplusplus)
target_link_libraries(Signal PUBLIC matplot)


# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(Signal
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
