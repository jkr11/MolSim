cmake_policy(SET CMP0079 NEW)

enable_testing()
include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}
)
file(GLOB_RECURSE MOLSIM_TEST
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

option(BUILD_TESTS "Build tests executable" OFF)

if (BUILD_TESTS)
    message("BUILD_TESTS is enabled")
    add_executable(MolSimTest ${MOLSIM_TEST}
            BoundaryConditionsTest.cpp
            XmlReaderTest.cpp
            XmlReaderTest.cpp
            TruncatedLennardJones.cpp
            IndexForceTest.cpp
            SingularGravityTest.cpp
            MixingTest.cpp)

    target_link_libraries(MolSimTest MolSimLib GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(MolSimTest)
endif ()
