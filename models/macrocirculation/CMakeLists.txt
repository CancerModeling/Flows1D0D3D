cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "Macrocirculation")
set(ProjectLib "Lib${ProjectName}")

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(src/macrocirculation SOURCES)
AUX_SOURCE_DIRECTORY(src/macrocirculation/communication SOURCES)

add_library(${ProjectLib} ${SOURCES})
target_link_libraries(${ProjectLib} PUBLIC Utils)
target_include_directories(${ProjectLib} INTERFACE ./src)
target_link_libraries(${ProjectLib} PRIVATE nlohmann_json::nlohmann_json
        nanoflann)

# ****************************************************************************
# Add apps
# ****************************************************************************
add_subdirectory(apps)

# ****************************************************************************
# Add unittests
# ****************************************************************************
option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")
if (ENABLE_UNIT_TESTS)
    enable_testing()

    add_subdirectory(unittests)

    find_package(Catch2 QUIET)
    # if cmake cannot find the package, we download it
    if (NOT Catch2_FOUND)
        # get catch2 for testing
        include(FetchContent)
        FetchContent_Declare(
                catch2_content
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG v2.13.1
        )
        FetchContent_MakeAvailable(catch2_content)
    endif()
endif()
