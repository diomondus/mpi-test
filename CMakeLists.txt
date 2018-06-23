cmake_minimum_required(VERSION 3.10)
project(mpi-test)
add_subdirectory(lab1)
add_subdirectory(lab2)
add_subdirectory(lab3)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

add_executable(mpi-test main.cpp)