#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#set(XERCESC_DIR /dss/dsshome1/lxc0B/di34miq/xercesc)
#include_directories(${XERCESC_DIR}/include)
#link_directories(${XERCESC_DIR}/lib64)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(ENABLE_OPENMP "Enable OpenMP support" OFF)
option(ENABLE_PARALLEL_STRATEGY_2 "Enable Parallelization Strategy 2" OFF)


if (ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif ()

##### Collect all cpp files except MolSim.cpp #####
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

list(REMOVE_ITEM MY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/MolSim.cpp")

add_library(MolSimLib
        ${MY_SRC}
        forces/HarmonicForce.cpp
        forces/HarmonicForce.h
        forces/TruncatedLennardJones.cpp
        forces/TruncatedLennardJones.h
        forces/IndexForce.cpp
        forces/IndexForce.h
)


target_link_libraries(MolSimLib
        PUBLIC
        spdlog::spdlog
)

if (ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(MolSimLib PUBLIC OpenMP::OpenMP_CXX)
endif ()

target_link_libraries(MolSimLib PRIVATE xerces-c)

target_compile_features(MolSimLib
        PUBLIC cxx_std_17
)

target_include_directories(MolSimLib PUBLIC
        ${libxsd_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/defs
        ${PROJECT_SOURCE_DIR}/forces
        ${PROJECT_SOURCE_DIR}/debug
        ${PROJECT_SOURCE_DIR}/calc
        ${PROJECT_SOURCE_DIR}/io
        ${PROJECT_SOURCE_DIR}/utils
)

add_executable(MolSim MolSim.cpp
        io/file/in/xml/input.cxx
        defs/Generators/ParticleGenerator.h
        defs/Generators/SpheroidGenerator.h
        defs/Generators/SpheroidGenerator.cpp
        forces/SingularGravity.h
        defs/Thermostat.cpp
        io/file/out/XmlWriter.h
        io/file/out/checkpoint-schema.hxx
        io/file/out/checkpoint-schema.cxx
        io/file/out/XmlWriter.cpp
        defs/Generators/MembraneGenerator.cpp
)

set(FILES_TO_DISABLE_WARNINGS
        io/file/in/xml/input.cxx
        io/file/in/xml/input.hxx
)

set_source_files_properties(${FILES_TO_DISABLE_WARNINGS}
        PROPERTIES
        COMPILE_OPTIONS "-w"
)

# find_package(OpenMP REQUIRED)


if (ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(MolSimLib PUBLIC OpenMP::OpenMP_CXX)
    if (ENABLE_PARALLEL_STRATEGY_2)
        target_compile_definitions(MolSimLib PUBLIC USE_PARALLEL_STRATEGY_2)
        target_compile_options(MolSimLib PRIVATE $<$<CXX_COMPILER_ID:GNU>:-DUSE_PARALLEL_STRATEGY_2>)
    endif ()
endif ()

target_link_libraries(MolSim PUBLIC MolSimLib)

target_compile_options(MolSimLib
        PUBLIC
        $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-nan-infinity-disabled>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -Wno-nan-infinity-disabled>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)



