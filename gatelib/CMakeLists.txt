add_library(Gate STATIC qhipster.c)
target_compile_definitions(Gate PUBLIC ACCELERATE_NEW_LAPACK)
target_include_directories(Gate
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
                  "/opt/homebrew/opt/libomp/include"
        PUBLIC "${CMAKE_SOURCE_DIR}/include"
                "${CMAKE_SOURCE_DIR}/statelib"
)

target_link_directories(Gate INTERFACE "/opt/homebrew/opt/libomp/lib")

target_link_libraries(Gate INTERFACE
        -lblas
        -llapack
        -lomp
        State
)

add_executable(gatelibTest test/testQhipster.c)

target_compile_options(gatelibTest PRIVATE
        -Xpreprocessor
        -fopenmp
        -O2)

add_test(Gate COMMAND ${gatelibTest})
add_dependencies(Gate gatelibTest)
add_custom_command(
        TARGET Gate
        COMMENT "Run tests on gatelib"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failures
)