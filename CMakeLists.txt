cmake_minimum_required(VERSION 3.27)
project(qlib C)
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

add_subdirectory(statelib)
add_subdirectory(gatelib)

set(QLIBSRC src/cplxutil.c src/linalg.c)
set(QLIBINC include)

set(CIRCLIBSRC circlib/src/diag.c circlib/src/circlib.c circlib/src/pauli.c)
set(CIRCLIBINC circlib/include)

set(MEASLIBSRC measlib/src/exact.c)
set(MEASLIBINC measlib/include)

add_subdirectory(extern/Unity)

# include double has to be called manually when building unity library
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)

add_compile_definitions(ACCELERATE_NEW_LAPACK)

set(circlibTest circlibTest)
add_executable(${circlibTest}
        circlib/test/testCirclib.c
        ${QLIBSRC}
        ${CIRCLIBSRC})

target_compile_options(${circlibTest} PRIVATE
        -Xpreprocessor
        -fopenmp)

target_include_directories(${circlibTest} PRIVATE
        test
        ${QLIBINC}
        ${CIRCLIBINC})

target_link_libraries(${circlibTest} PRIVATE
        -lblas
        -llapack
        -lomp
        State
        Gate
        unity
)

set(measlibTest measlibTest)
add_executable(${measlibTest}
        measlib/test/testExact.c
        ${QLIBSRC}
        ${CIRCLIBSRC}
        ${MEASLIBSRC})

target_compile_options(${measlibTest} PRIVATE
        -Xpreprocessor
        -fopenmp)

target_include_directories(${measlibTest} PRIVATE
        test
        ${QLIBINC}
        ${CIRCLIBINC}
        ${MEASLIBINC})

target_link_libraries(${measlibTest} PRIVATE
        -lblas
        -llapack
        -lomp
        State
        Gate
        unity
)


add_test(NAME ${circlibTest} COMMAND ${circlibTest})
add_test(NAME ${measlibTest} COMMAND ${measlibTest})

add_library(q SHARED
        ${QLIBSRC}
        ${STATELIBSRC}
        ${GATELIBSRC}
        ${CIRCLIBSRC}
        ${MEASLIBSRC})

target_compile_options(q PRIVATE
        -Xpreprocessor
        -fopenmp
        -O2)

target_include_directories(q PUBLIC
        ${QLIBINC}
        ${GATELIBINC}
        ${CIRCLIBINC}
        ${MEASLIBINC})

target_link_directories(q PUBLIC ${GATELIBLINK})

target_link_libraries(q
        -lBLAS
        -lLAPACK
        -lomp
        State
)

add_dependencies(q circlibTest)
add_dependencies(q measlibTest)

add_custom_command(
        TARGET q
        COMMENT "Run tests"
        DEPENDS
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failures
)
