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
            CLArgumentParserTest.cpp
            DirectSumContainerTest.cpp
            LennardJonesTest.cpp
            VerletIntegratorTest.cpp
            BoundaryConditionsTest.cpp
            PeriodicBoundaryTest.cpp
            PeriodicBoundary3DTest.cpp
            PeriodicAndReflectiveCornerTest.cpp
            testUtil.h
            XmlReaderTest.cpp
            ThermostatTest.cpp
            StatisticsTest.cpp
            WallTest.cpp
            CheckpointTest.cpp)

    target_link_libraries(MolSimTest MolSimLib GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(MolSimTest)
endif ()
