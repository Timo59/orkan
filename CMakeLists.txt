cmake_minimum_required(VERSION 3.27)
project(qlib C)

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

add_subdirectory(extern/Unity)

add_executable(gatelibTest
               gatelib/test/testQhipster.c
               src/cplxutil.c
               src/linalg.c
               gatelib/src/qhipster.c)

target_include_directories(gatelibTest PRIVATE
                           test
                           include
                           gatelib/include)

target_link_libraries(gatelibTest PRIVATE unity)


add_library(libq STATIC
            src/cplxutil.c
            src/linalg.c
            statelib/src/pstate.c
            gatelib/src/qhipster.c
            circlib/src/diag.c
            circlib/src/operator.c
            circlib/src/pauli.c
            measlib/src/exact.c)

target_include_directories(libq PRIVATE
                           include
                           statelib/include
                           gatelib/include
                           circlib/include
                           measlib/include)

