set(EXAMPLE_SOURCES render_pgm.c locality_demo.c)

foreach(src ${EXAMPLE_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE curves m)
    target_compile_options(${name} PRIVATE -O2)
    curves_configure_target(${name})
endforeach()
