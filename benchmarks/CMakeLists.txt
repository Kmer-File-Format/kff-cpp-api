
# add the executable
add_executable(reader_bench
    reader.cpp
)

target_link_libraries(reader_bench
    PRIVATE kff
)

add_executable(api_bench
    api_bench.cpp
)

target_link_libraries(api_bench
    PRIVATE kff
)
