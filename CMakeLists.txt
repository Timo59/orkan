cmake_minimum_required(VERSION 3.27)
project(Q C)
include(CTest)
set(CMAKE_C_STANDARD 17)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(STATUS "Building on macOS")
    add_definitions(-DMACOS)
    set(MACOS 1)
    set(LINUX 0)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Building on Linux")
    add_definitions(-DLINUX)
    set(MACOS 0)
    set(LINUX 1)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
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

add_subdirectory(extern/Unity)
add_subdirectory(internal)
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)                       # include double has to be called
                                                                                    # manually when building unity lib
add_library(BLAS_compiler_flags INTERFACE)                                          # Interface library specifying usage
target_link_directories(BLAS_compiler_flags                                         # of BLAS routines
        INTERFACE   $<${LINUX}:/usr/lib/x86_64-linux-gnu/openblas-pthread>
)
target_link_libraries(BLAS_compiler_flags INTERFACE
        $<${MACOS}:-lblas>
        $<${LINUX}:-lopenblas>
        $<${LINUX}:-lpthread>
        -llapack
)

add_library(OMP_compiler_flags INTERFACE)                                           # Interface library specifying usage
target_compile_options(OMP_compiler_flags INTERFACE                                 # requirements for dependents of
        $<${MACOS}:-Xpreprocessor>                                                  # the OpenMP API
        -fopenmp
)
target_include_directories(OMP_compiler_flags INTERFACE $<${MACOS}:/opt/homebrew/opt/libomp/include>)
target_link_libraries(OMP_compiler_flags INTERFACE $<${MACOS}:-lomp>)
target_link_directories(OMP_compiler_flags INTERFACE $<${MACOS}:/opt/homebrew/opt/libomp/lib>)

add_library(Q_compiler_flags INTERFACE)                                             # Interface library specifying usage
target_compile_options(Q_compiler_flags INTERFACE -O0 -g)                           # requirements for all Q dependents

add_library(Test_flags INTERFACE)                                                   # Interface library specifying usage
target_include_directories(Test_flags                                               # requirements for unit tests
        INTERFACE   ${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/internal/include
)
target_link_libraries(Test_flags INTERFACE
        Internal
        unity
)

add_subdirectory(state)
add_subdirectory(gate)
add_subdirectory(circ)
add_subdirectory(meas)
#add_executable(measlibTest meas/test/testExact.c)
#target_include_directories(measlibTest PRIVATE "${CMAKE_SOURCE_DIR}/include")
#target_link_libraries(measlibTest
#        PRIVATE
#            Test_flags
#            Math
#            State
#            Gate
#            Circ
#            Meas
#)

add_test(NAME testGate COMMAND gateTest)
add_test(NAME testCirc COMMAND circTest)
#add_test(NAME testMeas COMMAND measTest)

add_library(${PROJECT_NAME} SHARED
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

add_dependencies(Q
        gateTest
        circTest
#        measlibTest
)

add_custom_command(
        TARGET Q
        COMMENT "Run tests on qlib"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
)