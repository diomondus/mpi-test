cmake_minimum_required(VERSION 3.10)
project(mpi-test)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

add_executable(mpi-test main.cpp)