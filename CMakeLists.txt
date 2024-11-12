cmake_minimum_required(VERSION 3.27)
project(Q C)
include(CTest)
set(CMAKE_C_STANDARD 17)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}")
        endif()
    endif()
endif()

add_subdirectory(mathlib)
add_subdirectory(test)
add_subdirectory(extern/Unity)
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)                       # include double has to be called
                                                                                    # manually when building unity lib
add_library(BLAS_compiler_flags INTERFACE)                                          # Interface library specifying usage
                                                                                    # requirements for dependents of
target_compile_definitions(BLAS_compiler_flags INTERFACE ACCELERATE_NEW_LAPACK)     # BLAS specs
target_link_libraries(BLAS_compiler_flags INTERFACE
        -lblas
        -llapack
)

add_library(OMP_compiler_flags INTERFACE)                                           # Interface library specifying usage
                                                                                    # requirements for dependents of
target_compile_options(OMP_compiler_flags INTERFACE                                 # the OpenMP API
        -Xpreprocessor
        -fopenmp
)
target_include_directories(OMP_compiler_flags INTERFACE "/opt/homebrew/opt/libomp/include")
target_link_libraries(OMP_compiler_flags INTERFACE -lomp)
target_link_directories(OMP_compiler_flags INTERFACE "/opt/homebrew/opt/libomp/lib")

add_library(Q_compiler_flags INTERFACE)                                             # Interface library specifying usage
target_compile_options(Q_compiler_flags INTERFACE -O0)                              # requirements for all Q dependents

add_library(Test_flags INTERFACE)                                                   # Interface library specifying usage
target_link_libraries(Test_flags INTERFACE Test unity)                              # requirements for unit tests

add_subdirectory(statelib)
add_subdirectory(gatelib)

add_executable(gatelibTest gatelib/test/testQhipster.c)
target_include_directories(gatelibTest PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(gatelibTest
        PRIVATE Gate
                Math
                State
                Test_flags
)

add_subdirectory(circlib)
add_executable(circlibTest circlib/test/testCirclib.c)
target_include_directories(circlibTest PRIVATE "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(circlibTest
        PRIVATE Test_flags
                Math
                State
                Gate
                Circ
)

add_subdirectory(measlib)

add_executable(measlibTest measlib/test/testExact.c)

target_include_directories(measlibTest PRIVATE "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(measlibTest
        PRIVATE
            Test_flags
            Math
            State
            Gate
            Circ
            Meas
)

add_test(NAME gateTest COMMAND gatelibTest)
add_test(NAME circTest COMMAND circlibTest)
add_test(NAME measTest COMMAND measlibTest)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl")
add_library(${PROJECT_NAME} SHARED
        $<$<IN_LIST:${CIRCMETH},OLD;OMP;PROFILE>:$<TARGET_OBJECTS:Math>>
        $<$<IN_LIST:${MEASMETH},OLD;OMP;PROFILE>:$<TARGET_OBJECTS:Math>>
        $<TARGET_OBJECTS:State>
        $<TARGET_OBJECTS:Gate>
        $<TARGET_OBJECTS:Circ>
        $<TARGET_OBJECTS:Meas>
)
target_link_libraries(${PROJECT_NAME}
        PUBLIC  State
                Gate
                Circ
                Meas
                $<$<IN_LIST:${CIRCMETH},OLD;OMP;PROFILE>:Math>
                $<$<IN_LIST:${MEASMETH},OLD;OMP;PROFILE>:Math>
)

install(TARGETS Q
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
)
install(FILES       "include/qLib.h" "include/qTypes.h"
        DESTINATION include)

add_dependencies(Q gatelibTest circlibTest measlibTest)

add_custom_command(
        TARGET Q
        COMMENT "Run tests on gatelib"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failures
)