set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)

add_library(Circ OBJECT ${SOURCE_DIR}/pqc.c)
target_include_directories(Circ
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PUBLIC  ${INCLUDE_DIR}
)
target_link_libraries(Circ
        PRIVATE State
                BLAS_compiler_flags
                OMP_compiler_flags
                Q_compiler_flags
                Test_flags
)

add_executable(circTest ${TEST_DIR}/testPQC.c)
target_link_libraries(circTest
        PRIVATE State
                Gate
                Circ
                Test_flags
)

install(FILES ${INCLUDE_DIR}/pqc.h DESTINATION include/qlib
)
