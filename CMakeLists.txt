cmake_minimum_required(VERSION 4.0)
project(CppFM)

set(CMAKE_CXX_STANDARD 20)

add_executable(CppFM main.cpp
        Model.cpp
        Model.hpp)
