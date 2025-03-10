add_library(State OBJECT src/pstate.c)

target_include_directories(State
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        PUBLIC  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

install(FILES include/pstate.h DESTINATION include/qlib)