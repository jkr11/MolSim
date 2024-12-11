cmake_minimum_required(VERSION 3.10)

# Project name and version
project(Contest)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)
add_compile_options(-O3 -ffast-math -march=native -std=c++17 -g)

file(GLOB_RECURSE SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_executable(Contest ${SOURCES})

target_include_directories(Contest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(Contest
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
)