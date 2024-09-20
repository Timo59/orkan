add_library(Gate STATIC qhipster.c)

target_include_directories(Gate
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Gate
        PRIVATE Q_compiler_flags
                State
)