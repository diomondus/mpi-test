cmake_minimum_required(VERSION 3.10)

project(reduce)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

set(SOURCE_FILES reduce.cpp)

add_executable(reduce ${SOURCE_FILES})