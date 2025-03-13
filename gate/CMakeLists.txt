set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)

add_library(Gate OBJECT ${SOURCE_DIR}/qhipster.c)
target_include_directories(Gate                                 # Include directories for the object library
        PRIVATE ${CMAKE_SOURCE_DIR}/include                     # qTypes only populates the include directories for Gate
        PUBLIC  ${INCLUDE_DIR}                                  # Include files from gatelib are propagated
)
target_link_libraries(Gate                                      # Link libraries depending on the Cache string
        PRIVATE State
        PUBLIC  Q_compiler_flags
                BLAS_compiler_flags
)

add_executable(gateTest
        ${TEST_DIR}/testQhipster.c
)
target_link_libraries(gateTest
        PRIVATE Gate
                State
                Test_flags
)

install(FILES ${INCLUDE_DIR}/qhipster.h DESTINATION include/qlib)