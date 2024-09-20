add_library(State OBJECT pstate.c)

target_include_directories(State
        PRIVATE   "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

install(FILES statelib.h DESTINATION include)