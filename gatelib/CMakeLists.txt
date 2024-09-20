add_library(Gate OBJECT qhipster.c $<TARGET_OBJECTS:State>)

target_include_directories(Gate
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Gate
        PRIVATE Q_compiler_flags
        State
)

install(FILES gatelib.h DESTINATION include)