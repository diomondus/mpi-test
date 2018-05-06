cmake_minimum_required(VERSION 3.10)

project(latency)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

set(SOURCE_FILES latency.cpp)

add_executable(latency ${SOURCE_FILES})