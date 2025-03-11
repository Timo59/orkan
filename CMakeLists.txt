cmake_minimum_required(VERSION 3.27)
project(Q C)
include(CTest)
set(CMAKE_C_STANDARD 17)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(STATUS "Building on macOS")
    add_definitions(-DMACOS)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Building on Linux")
    add_definitions(-DLINUX)
else()
    message(FATAL_ERROR "Unknown system: ${CMAKE_SYSTEM_NAME}")
endif ()

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

set(CIRCMETH "BLAS" CACHE STRING "Choose circuit method of build: BLAS, BLAS_OMP, OLD, OMP, or PROFILE")
set(MEASMETH "BLAS_OMP" CACHE STRING "Choose measurement method of build: BLAS, BLAS_OMP, OLD, OMP, or PROFILE")
set_property(CACHE CIRCMETH PROPERTY STRINGS "BLAS" "BLAS_OMP" "OLD" "OMP" "OLD" "PROFILE")
set_property(CACHE MEASMETH PROPERTY STRINGS "BLAS" "BLAS_OMP" "OLD" "OMP" "OLD" "PROFILE")

add_subdirectory(extern/Unity)
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)                       # include double has to be called
                                                                                    # manually when building unity lib
add_library(BLAS_compiler_flags INTERFACE)                                          # Interface library specifying usage
target_link_libraries(BLAS_compiler_flags INTERFACE                                 # requirements for dependents of
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

set(TEST_SRC
        ${CMAKE_SOURCE_DIR}/internal/src/finDiff.c
        ${CMAKE_SOURCE_DIR}/internal/src/gateMat.c
        ${CMAKE_SOURCE_DIR}/internal/src/linalg.c
        ${CMAKE_SOURCE_DIR}/internal/src/pauliMat.c
        ${CMAKE_SOURCE_DIR}/internal/src/utils.c
)
add_library(Test_flags INTERFACE)                                                   # Interface library specifying usage
target_include_directories(Test_flags                                               # requirements for unit tests
        INTERFACE   ${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/internal/include
)
target_link_libraries(Test_flags INTERFACE unity)

add_subdirectory(state)
add_subdirectory(gate)

add_executable(gatelibTest
        gate/test/testQhipster.c
        ${TEST_SRC}
)
target_link_libraries(gatelibTest
        PRIVATE Gate
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
        INCLUDES DESTINATION include/qlib
)
install(FILES       "include/qLib.h"
                    "include/qTypes.h"
        DESTINATION include/qlib)

add_dependencies(Q gatelibTest circlibTest measlibTest)

add_custom_command(
        TARGET Q
        COMMENT "Run tests on gatelib"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
)