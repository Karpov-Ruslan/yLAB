cmake_minimum_required(VERSION 3.8)
project(HW3D)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O2")


add_executable(HW3D main.cpp geometry.cpp)

add_executable(Tests tests.cpp geometry.cpp)