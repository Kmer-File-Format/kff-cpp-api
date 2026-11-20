cmake_minimum_required(VERSION 3.16...3.31)

# The version of the project is the version of the library (its interface and its code). It is
# defined in kff_io.hpp (KFF_API_VERSION_MAJOR/MINOR), with the version of the KFF format written
# in the files (KFF_FORMAT_VERSION_MAJOR/MINOR, independent of the version of the library).
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/kff_io.hpp" KFF_API_VERSION_LINES
     REGEX "^#define KFF_API_VERSION_(MAJOR|MINOR) +[0-9]+ *$")
string(REGEX REPLACE ".*KFF_API_VERSION_MAJOR +([0-9]+).*" "\\1" KFF_API_MAJOR "${KFF_API_VERSION_LINES}")
string(REGEX REPLACE ".*KFF_API_VERSION_MINOR +([0-9]+).*" "\\1" KFF_API_MINOR "${KFF_API_VERSION_LINES}")
if (NOT KFF_API_MAJOR MATCHES "^[0-9]+$" OR NOT KFF_API_MINOR MATCHES "^[0-9]+$")
    message(FATAL_ERROR "Version of the library not found in kff_io.hpp")
endif()
project(KFF VERSION ${KFF_API_MAJOR}.${KFF_API_MINOR} LANGUAGES CXX)

# True when this directory is the main project (not added by add_subdirectory)
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(KFF_IS_TOP_LEVEL ON)
else()
    set(KFF_IS_TOP_LEVEL OFF)
endif()

option(KFF_BUILD_EXAMPLES "Build the kff_example program" ${KFF_IS_TOP_LEVEL})
option(KFF_BUILD_BENCHMARKS "Build the benchmark programs" ${KFF_IS_TOP_LEVEL})
option(KFF_INSTALL "Generate the install rules" ${KFF_IS_TOP_LEVEL})
option(KFF_BUILD_TESTS "Build the unit tests (run them with ctest)" ${KFF_IS_TOP_LEVEL})

# Optimised build by default (a parent project keeps control of its own type)
if (KFF_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

include(GNUInstallDirs)

add_library(kff kff_io.cpp)
add_library(kff::kff ALIAS kff)
set_target_properties(kff PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# The library is compiled in C++17; the public header only needs C++11.
target_compile_features(kff PRIVATE cxx_std_17 INTERFACE cxx_std_11)
set_target_properties(kff PROPERTIES CXX_EXTENSIONS OFF)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(kff PRIVATE -Wall -Wextra)
endif()

target_include_directories(kff
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if (KFF_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS kff EXPORT kffTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    install(FILES ${PROJECT_SOURCE_DIR}/kff_io.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # CMake package: find_package(kff) then target_link_libraries(... kff::kff)
    set(KFF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/kff)
    install(EXPORT kffTargets NAMESPACE kff:: DESTINATION ${KFF_CMAKE_DIR})
    configure_package_config_file(
        ${PROJECT_SOURCE_DIR}/cmake/kffConfig.cmake.in
        ${PROJECT_BINARY_DIR}/kffConfig.cmake
        INSTALL_DESTINATION ${KFF_CMAKE_DIR}
    )
    write_basic_package_version_file(
        ${PROJECT_BINARY_DIR}/kffConfigVersion.cmake
        COMPATIBILITY SameMajorVersion
    )
    install(FILES
        ${PROJECT_BINARY_DIR}/kffConfig.cmake
        ${PROJECT_BINARY_DIR}/kffConfigVersion.cmake
        DESTINATION ${KFF_CMAKE_DIR}
    )
endif()

if (KFF_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if (KFF_BUILD_BENCHMARKS)
    add_subdirectory("benchmarks")
endif()

if (KFF_BUILD_EXAMPLES)
    add_executable(kff_example main.cpp)
    target_link_libraries(kff_example PRIVATE kff)
endif()
