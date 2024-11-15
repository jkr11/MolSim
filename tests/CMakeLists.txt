enable_testing()
include_directories(
        ${PROJECT_SOURCE_DIR}/src
)
file(GLOB_RECURSE MOLSIM_TEST
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

option(BUILD_TESTS "Build tests executable" OFF)

if (BUILD_TESTS)
    add_executable(MolSimTest ${MOLSIM_TEST})


    target_link_libraries(MolSimTest MOlSimLib GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(MolSimTest)
endif ()
