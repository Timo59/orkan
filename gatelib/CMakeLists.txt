add_library(Gate STATIC qhipster.c)

target_include_directories(Gate
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
)

target_link_libraries(Gate
        PUBLIC  Q_compiler_flags
        PRIVATE State
)