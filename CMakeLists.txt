cmake_minimum_required(VERSION 3.14)
project(curves VERSION 1.1.0 LANGUAGES C)

option(CURVES_BUILD_TESTS        "Build test suite" ON)
option(CURVES_BUILD_BENCHMARKS   "Build benchmarks" OFF)
option(CURVES_BUILD_EXAMPLES     "Build examples"   ON)
option(CURVES_ENABLE_SANITIZERS  "ASan + UBSan on tests/examples/benchmarks" OFF)

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/curvesHelpers.cmake")

add_library(curves INTERFACE)
add_library(curves::curves ALIAS curves)
target_include_directories(curves INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(curves INTERFACE c_std_11)

if(CURVES_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(CURVES_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(CURVES_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

include(GNUInstallDirs)

install(FILES include/curves.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS curves EXPORT curvesTargets)
install(EXPORT curvesTargets
    FILE curvesTargets.cmake
    NAMESPACE curves::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curves
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/curvesConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/curvesConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curves
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/curvesConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/curvesConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/curvesConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/curves
)
