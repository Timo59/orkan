set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(State OBJECT ${SOURCE_DIR}/pstate.c)

target_include_directories(State
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PUBLIC  ${INCLUDE_DIR}
)
target_link_libraries(State
        PRIVATE BLAS_compiler_flags
        PUBLIC  $<${LINUX}:-lm>
)

install(FILES include/pstate.h DESTINATION include/qlib)