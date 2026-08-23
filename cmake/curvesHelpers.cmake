# Shared warning / sanitizer flags for tests, examples, and benchmarks.
# Not applied to the INTERFACE library — consumers choose their own flags.

function(curves_set_warnings target)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:C>:-Wall;-Wextra;-Wpedantic;-Wshadow>
            $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wshadow>
        )
    elseif(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    endif()
endfunction()

function(curves_set_sanitizers target)
    if(NOT CURVES_ENABLE_SANITIZERS)
        return()
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -fno-sanitize-recover=address,undefined
        )
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined
        )
    endif()
endfunction()

function(curves_configure_target target)
    curves_set_warnings(${target})
    curves_set_sanitizers(${target})
endfunction()
