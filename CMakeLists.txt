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
# include double has to be called manually when building unity library
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)

# create an interface library to propagate linking libraries for all tests
add_library(Test_flags INTERFACE)
target_link_libraries(Test_flags INTERFACE Math Test unity)

# create an interface library to propagate compiler flags as well as include and link directories and link libraries to
# all subdirectories using LAPACK and/or omp
add_library(Q_compiler_flags INTERFACE)
target_compile_definitions(Q_compiler_flags INTERFACE ACCELERATE_NEW_LAPACK)
target_compile_options(Q_compiler_flags INTERFACE
        -Xpreprocessor
        -fopenmp
        -O2
)

target_include_directories(Q_compiler_flags INTERFACE "/opt/homebrew/opt/libomp/include")
target_link_libraries(Q_compiler_flags INTERFACE
        -lblas
        -llapack
        -lomp
)

target_link_directories(Q_compiler_flags INTERFACE "/opt/homebrew/opt/libomp/lib")

add_subdirectory(statelib)
add_subdirectory(gatelib)

add_executable(gatelibTest gatelib/test/testQhipster.c)

target_include_directories(gatelibTest PRIVATE "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(gatelibTest
        PRIVATE Test_flags
        State
        Gate
)

add_test(gateTest COMMAND gatelibTest)

add_library(Q SHARED qlib.c)
target_link_libraries(Q PRIVATE Math)

add_dependencies(Q gatelibTest)

add_custom_command(
        TARGET Q
        COMMENT "Run tests on gatelib"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failures
)