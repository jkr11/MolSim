cmake_minimum_required(VERSION 3.10)

project(PSEMolDyn_GroupA VERSION 0.0.3 LANGUAGES CXX)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-O3 -ffast-math -march=native -std=c++17 -Wno-deprecated-declarations -fopenmp)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
include(googletest)
include(doxygen)
include(spdlog)
include(xsd)

# check for memory leaks with asan (run time) and clang-tidy (compile time)
if (CMAKE_BUILD_TYPE STREQUAL "asan")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address -g")
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=clang-analyzer-*;-extra-arg=-std=c++17")
    execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean)
    message("Check mode is enabled (sanitize address, clang-tidy), clean build")
endif ()

# check for memory leaks with asan (run time) and clang-tidy (compile time), for github
if (CMAKE_BUILD_TYPE STREQUAL "asan-quiet")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address -g")
    set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=clang-analyzer-*;-std=c++17;--quiet;-extra-arg=-std=c++17")
    execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean)
    message("Check mode is enabled (sanitize address, clang-tidy), clean build, quiet")
endif ()

# Enable DEBUG macro
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-ggdb3)
    add_compile_definitions(DEBUG)
    message("Debug mode is enabled")
endif ()

option(ENABLE_BENCHMARK "Enable benchmark mode" ON)

if (ENABLE_BENCHMARK)
    message("Benchmark mode is enabled")
    add_compile_definitions(BENCHMARK)
endif ()

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_subdirectory(src)
add_subdirectory(tests)
