cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

project(metodyObliczeniowe)

include_directories(include)

add_executable(main src/main.cpp src/calerf.cpp src/func.cpp)