set(TEST_SOURCES
    test_curves.c
    test_wide_keys.c
    test_prefix.c
    test_range.c
    test_quant.c
    test_simd.c
    test_research.c
    test_dispatch.c
)

foreach(src ${TEST_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE curves m)
    target_compile_options(${name} PRIVATE -O2)
    curves_configure_target(${name})
    add_test(NAME ${name} COMMAND ${name})
endforeach()

enable_language(CXX)
add_executable(test_cxx test_cxx.cpp)
target_link_libraries(test_cxx PRIVATE curves)
target_compile_features(test_cxx PRIVATE cxx_std_11)
target_compile_options(test_cxx PRIVATE -O2)
curves_configure_target(test_cxx)
add_test(NAME test_cxx COMMAND test_cxx)
