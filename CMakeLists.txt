cmake_minimum_required(VERSION 3.4)
project(MiniMD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)

set(SOURCE_FILES main.cpp Simulation.cpp Simulation.hpp helper.cpp helper.hpp)
add_executable(MiniMD ${SOURCE_FILES})