add_executable(kff_tests
    main.cpp
    test_header.cpp
    test_variables.cpp
    test_raw.cpp
    test_minimizer.cpp
    test_index.cpp
    test_reader.cpp
    test_file.cpp
)
target_link_libraries(kff_tests PRIVATE kff)
target_compile_features(kff_tests PRIVATE cxx_std_17)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(kff_tests PRIVATE -Wall -Wextra)
endif()

add_test(NAME kff_tests COMMAND kff_tests)
