cmake_policy(SET CMP0135 NEW)
include(FetchContent)
FetchContent_Declare(
        libxsd
        URL
        ${PROJECT_SOURCE_DIR}/libs/libxsd-4.0.0.zip
        #URL_HASH MD5=4d4832a6bc30cd96e1e119c47d9ca3e0
)

FETCHCONTENT_MAKEAVAILABLE(libxsd)