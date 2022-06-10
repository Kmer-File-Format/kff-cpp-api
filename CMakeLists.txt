cmake_minimum_required(VERSION 3.19)

# set the project name
project(KFF)
project(KFF VERSION 1.0)
configure_file(kff_io.hpp.in kff_io.hpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_BUILD_TYPE Debug)

include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "-Wall -Ofast -g")

#include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(kff kff_io.cpp)

target_include_directories(kff
    PUBLIC
    .
    ${CMAKE_CURRENT_SOURCE_DIR}/
)


# add the executable
add_executable(kff_example
    main.cpp
)

target_link_libraries(kff_example
    kff
)

