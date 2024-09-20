add_library(Circ STATIC diag.c pauli.c circlib.c)

target_include_directories(Circ
        PRIVATE   "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Circ PRIVATE
        Q_compiler_flags
        Math
        State
        Gate
)