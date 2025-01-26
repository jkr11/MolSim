include_directories(${CMAKE_CURRENT_SOURCE_DIR})

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
        xerces-c
        spdlog::spdlog
)

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

target_link_libraries(MolSim PUBLIC MolSimLib)

target_compile_options(MolSimLib
        PUBLIC
        $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-nan-infinity-disabled>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -Wno-nan-infinity-disabled>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

