cmake_minimum_required(VERSION 3.21)
project(QuantumGravity)

set(CMAKE_CXX_STANDARD 14)

add_executable(QuantumGravity main.cpp
        src/LatticeFieldTheory.cpp
        src/LatticeFieldTheory.h
        src/Point.cpp src/Point.h tests.cpp tests.h)
