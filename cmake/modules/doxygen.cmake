# option to disable doxygen
option(ENABLE_DOXYGEN "Enable doxygen" ON)

if (NOT ENABLE_DOXYGEN)
    message("Doxygen is disabled")
    return()
endif ()

# find doxygen install
find_package(Doxygen)

if (NOT DOXYGEN_FOUND)
    message("No doxygen install found")
    return()
endif ()

# custom doc_doxygen target
add_custom_target(doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        VERBATIM
)

set_target_properties(doc_doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)