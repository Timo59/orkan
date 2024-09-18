add_library(State STATIC pstate.c)
target_include_directories(State INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
)