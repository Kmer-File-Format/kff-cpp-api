cmake_minimum_required(VERSION 3.16...3.31)

# set the project name
project(KFF VERSION 1.0 LANGUAGES CXX)

# True when this directory is the main project (not added by add_subdirectory)
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(KFF_IS_TOP_LEVEL ON)
else()
    set(KFF_IS_TOP_LEVEL OFF)
endif()

# Optimised build by default (a parent project keeps control of its own type)
if (KFF_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# The generated header is placed in this project's build directory, so that
# the library also works when it is included with add_subdirectory().
configure_file(
    ${PROJECT_SOURCE_DIR}/kff_io.hpp.in
    ${PROJECT_BINARY_DIR}/kff_io.hpp
)

add_library(kff kff_io.cpp)
add_library(kff::kff ALIAS kff)

# The library is compiled in C++17; the public header only needs C++11.
target_compile_features(kff PRIVATE cxx_std_17 INTERFACE cxx_std_11)
set_target_properties(kff PROPERTIES CXX_EXTENSIONS OFF)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(kff PRIVATE -Wall -Wextra)
endif()

target_include_directories(kff
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)

add_subdirectory("benchmarks")

# add the executable
add_executable(kff_example
    main.cpp
)

target_link_libraries(kff_example
    kff
)
