cmake_policy(SET CMP0135 NEW)
include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.14.1
)

FETCHCONTENT_MAKEAVAILABLE(spdlog)