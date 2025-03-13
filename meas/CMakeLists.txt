set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Meas OBJECT ${SOURCE_DIR}/exact.c)
target_include_directories(Meas
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PUBLIC  ${INCLUDE_DIR}
)
target_link_libraries(Meas
        PRIVATE State
                Circ
                BLAS_compiler_flags
                OMP_compiler_flags
                Q_compiler_flags
)

install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/exact.h
        DESTINATION include/qlib)