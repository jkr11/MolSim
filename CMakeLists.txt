cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn_GroupA VERSION 0.0.1)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

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

# debug
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG)
    message("Debug mode is enabled")
endif ()

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

##### Collect all cpp files except MolSim.cpp #####
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

list(REMOVE_ITEM MY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/MolSim.cpp")


##### fetch spdlog and googletest #####
include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.14.1
)

FETCHCONTENT_MAKEAVAILABLE(spdlog)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.1
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)


##### Build MolSimLib #####
add_library(MolSimLib
        ${MY_SRC}
)

target_link_libraries(MolSimLib
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
        spdlog::spdlog
)

target_include_directories(MolSimLib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_features(MolSimLib
        PUBLIC cxx_std_17
)

target_compile_options(MolSimLib
        PUBLIC
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")


##### Build MolSim.exe #####
add_executable(MolSim
        src/MolSim.cpp
)

target_link_libraries(MolSim
        PRIVATE MolSimLib
)

##### Build MolSimTest.exe #####
enable_testing()

add_executable(MolSimTest
    tests/ParticleContainerTest.cpp
    tests/VerletIntegratorTest.cpp
    tests/LennardJonesTest.cpp
    tests/CLArgumentParserTest.cpp
)

target_link_libraries(MolSimTest
        PRIVATE MolSimLib
        PUBLIC GTest::gtest_main
)

gtest_discover_tests(MolSimTest)


##### Other #####
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)