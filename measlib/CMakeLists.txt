add_library(Meas STATIC exact.c)

target_include_directories(Meas
        PRIVATE   "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Meas PRIVATE
        Q_compiler_flags
        Math
        State
        Gate
        Circ
)