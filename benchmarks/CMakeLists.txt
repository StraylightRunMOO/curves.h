add_executable(bench_main bench_main.c)
target_link_libraries(bench_main PRIVATE curves m)
target_compile_options(bench_main PRIVATE -O3 -march=native)
curves_configure_target(bench_main)
